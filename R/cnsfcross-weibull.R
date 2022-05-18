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
# Convolution: weibull - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf weibull-normal distribution
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
# same sigma_u

## logit specification class membership
ccnsfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
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
  Pi2 <- numeric(N)
  if (k < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv2[i]/2)))
    }
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cauchit specification class membership
ccnsfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  Pi2 <- numeric(N)
  if (k < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv2[i]/2)))
    }
    Probc1 <- 1/pi * atan(Wz) + 1/2
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## probit specification class membership
ccnsfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
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
  Pi2 <- numeric(N)
  if (k < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv2[i]/2)))
    }
    Probc1 <- pnorm(Wz)
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cloglog specification class membership
ccnsfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  Pi2 <- numeric(N)
  if (k < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur)/exp(Wv2[i]/2)))
    }
    Probc1 <- 1 - exp(-exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# different sigma_u

## logit specification class membership
cmcesfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur1)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur2)/exp(Wv2[i]/2)))
    }
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cauchit specification class membership
cmcesfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur1)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur2)/exp(Wv2[i]/2)))
    }
    Probc1 <- 1/pi * atan(Wz) + 1/2
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## probit specification class membership
cmcesfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur1)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur2)/exp(Wv2[i]/2)))
    }
    Probc1 <- pnorm(Wz)
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cloglog specification class membership
cmcesfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(NA)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        S * ur1)/exp(Wv1[i]/2)))
      Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        S * ur2)/exp(Wv2[i]/2)))
    }
    Probc1 <- 1 - exp(-exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# starting value for the log-likelihood ----------
#' starting values for cnsf weibull-normal distribution
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
# same sigma_u
cstcnsfweibullnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
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
      1]), Yvar = Yvar, Xvar = Xvar, S = S, N = N, FiMat = FiMat,
    wHvar = wHvar)
  Esti <- initWeibull$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), "k", paste0("CN_", colnames(Zvar)))
  names(initWeibull$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "k")
  return(list(StartVal = StartVal, initWeibull = initWeibull))
}

# different sigma_u
cstmcesfweibullnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + weibull - normal distributions...\n")
  initWeibull <- maxLik(logLik = cweibullnormlike, start = cstweibullnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradweibullnormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = printInfo, reltol = tol), nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]),
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat)
  Esti <- initWeibull$estimate
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar + 3], 1.05 *
    Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), "k",
    "k", paste0("MCE_", colnames(Zvar)))
  names(initWeibull$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "k")
  return(list(StartVal = StartVal, initWeibull = initWeibull))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf weibull-normal distribution
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
# same sigma_u

## logit specification class membership
cgradcnsfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
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
  ewu_h <- exp(Wu/2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewu_h, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"))
  F2_2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewv2_h, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1, MARGIN = 1, STATS = 1/ewv1_h^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1, MARGIN = 1, STATS = 1/ewv2_h^3,
    FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewv2_h,
    FUN = "*"), 1, sum)
  sigx1 <- (prC * sDiv2 + ewz * sDiv1/wzdeno)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * F1, MARGIN = 1,
    STATS = ewu_h/ewv1_h^3, FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * F1, MARGIN = 1,
    STATS = ewu_h/ewv2_h^3, FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1^2), MARGIN = 1, STATS = 1/ewv1_h^2,
    FUN = "*") - 0.5 * F2_1, MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1^2), MARGIN = 1, STATS = 1/ewv2_h^2,
    FUN = "*") - 0.5 * F2_2, MARGIN = 1, STATS = 1/ewv2_h,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewu_h/(k^2 * ewv1_h^3),
    FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewu_h/(k^2 * ewv2_h^3),
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz/(sigx1 *
      wzdeno)
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gu <- gu1 + gu2
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2, 1, sum) * prC/sigx1 +
    apply(F6_1, 1, sum) * ewz/(sigx1 * wzdeno), sweep(Zvar,
    MARGIN = 1, STATS = prC * ewz * (sDiv1 - sDiv2)/(sigx1 *
      wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradcnsfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F2_2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr2, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1, MARGIN = 1, STATS = 1/ewvsr1^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1, MARGIN = 1, STATS = 1/ewvsr2^3,
    FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*"), 1, sum)
  sigx1 <- (ewz2 * sDiv2 + ewz1 * sDiv1)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr1^3, FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr2^3, FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1^2), MARGIN = 1, STATS = 1/ewvsr1^2,
    FUN = "*") - 0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1^2), MARGIN = 1, STATS = 1/ewvsr2^2,
    FUN = "*") - 0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewusr/(k^2 * ewvsr1^3),
    FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewusr/(k^2 * ewvsr2^3),
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz1/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gu <- gu1 + gu2
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2, 1, sum) * ewz2/sigx1 +
    apply(F6_1, 1, sum) * ewz1/sigx1, sweep(Zvar, MARGIN = 1,
    STATS = (sDiv1 - sDiv2)/(pi * sigx1 * ((Wz)^2 + 1)),
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradcnsfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
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
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F2_2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr2, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1, MARGIN = 1, STATS = 1/ewvsr1^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1, MARGIN = 1, STATS = 1/ewvsr2^3,
    FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*"), 1, sum)
  sigx1 <- ((1 - pwZ) * sDiv2 + pwZ * sDiv1)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr1^3, FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr2^3, FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1^2), MARGIN = 1, STATS = 1/ewvsr1^2,
    FUN = "*") - 0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1^2), MARGIN = 1, STATS = 1/ewvsr2^2,
    FUN = "*") - 0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewusr/(k^2 * ewvsr1^3),
    FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewusr/(k^2 * ewvsr2^3),
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * pwZ/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gu <- gu1 + gu2
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2, 1, sum) * (1 -
    pwZ)/sigx1 + apply(F6_1, 1, sum) * pwZ/sigx1, sweep(Zvar,
    MARGIN = 1, STATS = dwZ * (sDiv1 - sDiv2)/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradcnsfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  F1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k)), MARGIN = 1,
    STATS = ewusr, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F2_2 <- dnorm(sweep(F1, MARGIN = 1, STATS = 1/ewvsr2, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1, MARGIN = 1, STATS = 1/ewvsr1^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1, MARGIN = 1, STATS = 1/ewvsr2^3,
    FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*"), 1, sum)
  sigx1 <- ((1 - prZ) * sDiv1 + prZ * sDiv2)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr1^3, FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * F1, MARGIN = 1,
    STATS = ewusr/ewvsr2^3, FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1^2), MARGIN = 1, STATS = 1/ewvsr1^2,
    FUN = "*") - 0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1^2), MARGIN = 1, STATS = 1/ewvsr2^2,
    FUN = "*") - 0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k) * F2_1 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewusr/(k^2 * ewvsr1^3),
    FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k) * F2_2 * log(-log(1 -
    FiMat)) * F1, MARGIN = 1, STATS = S * ewusr/(k^2 * ewvsr2^3),
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gu <- gu1 + gu2
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2, 1, sum) * prZ/sigx1 +
    apply(F6_1, 1, sum) * (1 - prZ)/sigx1, sweep(Zvar, MARGIN = 1,
    STATS = prZ * ewz * (sDiv1 - sDiv2)/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u

## logit specification class membership
cgradmcesfweibullnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  F1_1 <- sweep(sweep(S * (-log(1 - FiMat))^(1/k1), MARGIN = 1,
    STATS = ewu1_h, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F1_2 <- sweep(sweep(S * (-log(1 - FiMat))^(1/k2), MARGIN = 1,
    STATS = ewu2_h, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewv2_h, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewv1_h^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewv2_h^3,
    FUN = "*")
  F4_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"))
  F4_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewv2_h, FUN = "*"))
  sDiv1 <- apply(sweep(F4_1, MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F4_2, MARGIN = 1, STATS = 1/ewv2_h,
    FUN = "*"), 1, sum)
  sigx1 <- (prC * sDiv2 + ewz * sDiv1/wzdeno)
  F5_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1,
    STATS = ewu1_h/ewv1_h^3, FUN = "*")
  F5_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1,
    STATS = ewu2_h/ewv2_h^3, FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 -
    FiMat)) * F1_1, MARGIN = 1, STATS = S * ewu1_h/(k1^2 *
    ewv1_h^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 -
    FiMat)) * F1_2, MARGIN = 1, STATS = S * ewu2_h/(k2^2 *
    ewv2_h^3), FUN = "*")
  F7_1 <- sweep(sweep(F2_1 * F1_1^2, MARGIN = 1, STATS = 0.5/ewv1_h^2,
    FUN = "*") - 0.5 * F4_1, MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  F7_2 <- sweep(sweep(F2_2 * F1_2^2, MARGIN = 1, STATS = 0.5/ewv2_h^2,
    FUN = "*") - 0.5 * F4_2, MARGIN = 1, STATS = 1/ewv2_h,
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
  }
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F5_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz/(sigx1 *
      wzdeno)
  }
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu2[, k] <- apply(sweep(-(0.5 * (S * F5_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
  }
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, apply(F6_1, 1, sum) *
    ewz/(sigx1 * wzdeno), apply(F6_2, 1, sum) * prC/sigx1,
    sweep(Zvar, MARGIN = 1, STATS = prC * ewz * (sDiv1 -
      sDiv2)/(sigx1 * wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmcesfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  F1_1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k1)), MARGIN = 1,
    STATS = ewusr1, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F1_2 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k2)), MARGIN = 1,
    STATS = ewusr2, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewvsr2, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewvsr1^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewvsr2^3,
    FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*"), 1, sum)
  sigx1 <- (ewz2 * sDiv2 + ewz1 * sDiv1)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1,
    STATS = ewusr1/ewvsr1^3, FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1,
    STATS = ewusr2/ewvsr2^3, FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1_1^2), MARGIN = 1, STATS = 1/ewvsr1^2,
    FUN = "*") - 0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1_2^2), MARGIN = 1, STATS = 1/ewvsr2^2,
    FUN = "*") - 0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 -
    FiMat)) * F1_1, MARGIN = 1, STATS = S * ewusr1/(k1^2 *
    ewvsr1^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 -
    FiMat)) * F1_2, MARGIN = 1, STATS = S * ewusr2/(k2^2 *
    ewvsr2^3), FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz1/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, apply(F6_1, 1, sum) *
    ewz1/sigx1, apply(F6_2, 1, sum) * ewz2/sigx1, sweep(Zvar,
    MARGIN = 1, STATS = (sDiv1 - sDiv2)/(pi * sigx1 * ((Wz)^2 +
      1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmcesfweibullnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  F1_1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k1)), MARGIN = 1,
    STATS = ewusr1, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F1_2 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k2)), MARGIN = 1,
    STATS = ewusr2, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewvsr2, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewvsr1^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewvsr2^3,
    FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*"), 1, sum)
  sigx1 <- ((1 - pwZ) * sDiv2 + pwZ * sDiv1)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1,
    STATS = ewusr1/ewvsr1^3, FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1,
    STATS = ewusr2/ewvsr2^3, FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1_1^2), MARGIN = 1, STATS = 1/ewvsr1^2,
    FUN = "*") - 0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1_2^2), MARGIN = 1, STATS = 1/ewvsr2^2,
    FUN = "*") - 0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 -
    FiMat)) * F1_1, MARGIN = 1, STATS = S * ewusr1/(k1^2 *
    ewvsr1^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 -
    FiMat)) * F1_2, MARGIN = 1, STATS = S * ewusr2/(k2^2 *
    ewvsr2^3), FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * pwZ/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, apply(F6_1, 1, sum) *
    pwZ/sigx1, apply(F6_2, 1, sum) * (1 - pwZ)/sigx1, sweep(Zvar,
    MARGIN = 1, STATS = dwZ * (sDiv1 - sDiv2)/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmcesfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  F1_1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k1)), MARGIN = 1,
    STATS = ewusr1, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F1_2 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k2)), MARGIN = 1,
    STATS = ewusr2, FUN = "*"), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewvsr1, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewvsr2, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewvsr1^3,
    FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewvsr2^3,
    FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*"), 1, sum)
  sigx1 <- ((1 - prZ) * sDiv1 + prZ * sDiv2)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1,
    STATS = ewusr1/ewvsr1^3, FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1,
    STATS = ewusr2/ewvsr2^3, FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1_1^2), MARGIN = 1, STATS = 1/ewvsr1^2,
    FUN = "*") - 0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr1,
    FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1_2^2), MARGIN = 1, STATS = 1/ewvsr2^2,
    FUN = "*") - 0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr2,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 -
    FiMat)) * F1_1, MARGIN = 1, STATS = S * ewusr1/(k1^2 *
    ewvsr1^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 -
    FiMat)) * F1_2, MARGIN = 1, STATS = S * ewusr2/(k2^2 *
    ewvsr2^3), FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, apply(F6_1, 1, sum) *
    (1 - prZ)/sigx1, apply(F6_2, 1, sum) * prZ/sigx1, sweep(Zvar,
    MARGIN = 1, STATS = prZ * ewz * (sDiv1 - sDiv2)/sigx1,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf weibull-normal distribution
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
# same sigma_u

## logit specification class membership
cnsfweibullnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfweibullnormlike_logit(startVal, nXvar = nXvar,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfweibullnormlike_logit,
      grad = cgradcnsfweibullnormlike_logit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(jacobian(function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), parm), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfweibullnormlike_logit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) -colSums(cgradcnsfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$solution)
  }
  mleObj$logL_OBS <- ccnsfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cauchit specification class membership
cnsfweibullnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfweibullnormlike_cauchit(startVal,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfweibullnormlike_cauchit,
      grad = cgradcnsfweibullnormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradcnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfweibullnormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- ccnsfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## probit specification class membership
cnsfweibullnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfweibullnormlike_probit(startVal,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfweibullnormlike_probit,
      grad = cgradcnsfweibullnormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradcnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfweibullnormlike_probit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- ccnsfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cloglog specification class membership
cnsfweibullnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfweibullnormlike_cloglog(startVal,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfweibullnormlike_cloglog,
      grad = cgradcnsfweibullnormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradcnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(ccnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfweibullnormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradcnsfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- ccnsfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

# different sigma_u

## logit specification class membership
mcesfweibullnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfweibullnormlike_logit(startVal,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfweibullnormlike_logit,
      grad = cgradmcesfweibullnormlike_logit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(jacobian(function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), parm), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfweibullnormlike_logit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) -colSums(cgradmcesfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$solution)
  }
  mleObj$logL_OBS <- cmcesfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfweibullnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cauchit specification class membership
mcesfweibullnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfweibullnormlike_cauchit(startVal,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfweibullnormlike_cauchit,
      grad = cgradmcesfweibullnormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmcesfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfweibullnormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmcesfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfweibullnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## probit specification class membership
mcesfweibullnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfweibullnormlike_probit(startVal,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfweibullnormlike_probit,
      grad = cgradmcesfweibullnormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmcesfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfweibullnormlike_probit(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmcesfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfweibullnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cloglog specification class membership
mcesfweibullnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfweibullnormlike_cloglog(startVal,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfweibullnormlike_cloglog,
      grad = cgradmcesfweibullnormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmcesfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmcesfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfweibullnormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmcesfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmcesfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfweibullnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}


# Conditional efficiencies estimation ----------
#' efficiencies for cnsf weibull-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfweibullnormeff_logit <- function(object, level) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
ccnsfweibullnormeff_cauchit <- function(object, level) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
ccnsfweibullnormeff_probit <- function(object, level) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
ccnsfweibullnormeff_cloglog <- function(object, level) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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

# different sigma_u

## logit specification class membership
cmcesfweibullnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k1, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur1)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur2)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmcesfweibullnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k1, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur1)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur2)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmcesfweibullnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k1, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur1)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur2)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmcesfweibullnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  u_c2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
      sigmaV = exp(Wv1[i]/2), k = k1, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
      sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
  }
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- numeric(object$Nobs)
    teBC_reciprocal_c1 <- numeric(object$Nobs)
    teBC_c2 <- numeric(object$Nobs)
    teBC_reciprocal_c2 <- numeric(object$Nobs)
    for (i in seq_along(1:object$Nobs)) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur1)/exp(Wv1[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv1[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv1[i]/2),
        k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i,
        ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur2)/exp(Wv2[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv2[i]/2), k = k2, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv2[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
#' marginal impact on efficiencies for cnsf weibull-normal distribution
#' @param object object of class sfacross
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfmargweibullnorm_Eu_logit <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargweibullnorm_Vu_logit <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
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
ccnsfmargweibullnorm_Eu_cauchit <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargweibullnorm_Vu_cauchit <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
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
ccnsfmargweibullnorm_Eu_probit <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargweibullnorm_Vu_probit <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
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
ccnsfmargweibullnorm_Eu_cloglog <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargweibullnorm_Vu_cloglog <- function(object) {
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
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2)))
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# different sigma_u

## logit specification class membership
cmcesfmargweibullnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2) * gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2) * gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmargweibullnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2),
      ncol = 1))
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
cmcesfmargweibullnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2) * gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2) * gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmargweibullnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2),
      ncol = 1))
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
cmcesfmargweibullnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2) * gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2) * gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmargweibullnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2),
      ncol = 1))
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
cmcesfmargweibullnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2) * gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2) * gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmargweibullnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur1)/exp(Wv1[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur2)/exp(Wv2[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2),
      ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2),
      ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}
