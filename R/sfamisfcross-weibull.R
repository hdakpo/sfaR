################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Multi-modal inefficiency (MISF)                                       #
# Link functions: - logit exp(theta * Z)/(1 + exp(theta * Z))                  #
#                 - cauchit 1/pi * atan(theta * Z) + 1/2                       #
#                 - probit pnorm(theta * Z)                                    #
#                 - cloglog 1 - exp(-exp(theta * Z))                           #
# Convolution: weibull - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for misf weibull-normal distribution
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
# logit specification class membership
cmisfweibullnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(-.Machine$double.xmax)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur1)/exp(Wv[i]/2)))
      Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur2)/exp(Wv[i]/2)))
    }
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar *
      log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# cauchit specification class membership
cmisfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(-.Machine$double.xmax)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur1)/exp(Wv[i]/2)))
      Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur2)/exp(Wv[i]/2)))
    }
    Probc1 <- 1/pi * atan(Wz) + 1/2
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar *
      log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# probit specification class membership
cmisfweibullnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(-.Machine$double.xmax)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur1)/exp(Wv[i]/2)))
      Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur2)/exp(Wv[i]/2)))
    }
    Probc1 <- pnorm(Wz)
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar *
      log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# cloglog specification class membership
cmisfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  if (k1 < 0 || k2 < 0) {
    return(-.Machine$double.xmax)
  } else {
    for (i in 1:N) {
      ur1 <- exp(Wu1[i]/2) * (-log(1 - FiMat[i, ]))^(1/k1)
      ur2 <- exp(Wu2[i]/2) * (-log(1 - FiMat[i, ]))^(1/k2)
      Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur1)/exp(Wv[i]/2)))
      Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S * ur2)/exp(Wv[i]/2)))
    }
    Probc1 <- 1 - exp(-exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar *
      log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# starting value for the log-likelihood ----------
#' starting values for misf weibull-normal distribution
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
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
cstmisfweibullnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat, whichStart, initIter, initAlg,
  printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- cstweibullnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initWeibull <- NULL
  } else {
    cat("Initialization: SFA + weibull - normal distributions...\n")
    initWeibull <- maxLik::maxLik(logLik = cweibullnormlike, start = cstweibullnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradweibullnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)
    Esti <- initWeibull$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1),
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar +
      3], 1.05 * Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "k", "k", paste0("MISF_",
    colnames(Zvar)))
  return(list(StartVal = StartVal, initWeibull = initWeibull))
}

# Gradient of the likelihood function ----------
#' gradient for misf weibull-normal distribution
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
## logit specification class membership
cgradmisfweibullnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv_h <- exp(Wv/2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  F1_1 <- sweep(sweep(S * (-log(1 - FiMat))^(1/k1), MARGIN = 1, STATS = ewu1_h,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F1_2 <- sweep(sweep(S * (-log(1 - FiMat))^(1/k2), MARGIN = 1, STATS = ewu2_h,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F4_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  F4_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  sDiv1 <- apply(sweep(F4_1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F4_2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"), 1, sum)
  sigx1 <- (prC * sDiv2 + ewz * sDiv1/wzdeno)
  F5_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1, STATS = ewu1_h/ewv_h^3,
    FUN = "*")
  F5_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1, STATS = ewu2_h/ewv_h^3,
    FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 - FiMat)) * F1_1,
    MARGIN = 1, STATS = S * ewu1_h/(k1^2 * ewv_h^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 - FiMat)) * F1_2,
    MARGIN = 1, STATS = S * ewu2_h/(k2^2 * ewv_h^3), FUN = "*")
  F7_1 <- sweep(sweep(F2_1 * F1_1^2, MARGIN = 1, STATS = 0.5/ewv_h^2, FUN = "*") -
    0.5 * F4_1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  F7_2 <- sweep(sweep(F2_2 * F1_2^2, MARGIN = 1, STATS = 0.5/ewv_h^2, FUN = "*") -
    0.5 * F4_2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * ewz/(sigx1 * wzdeno)
  }
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * prC/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F5_1)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
  }
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu2[, k] <- apply(sweep(-(0.5 * (S * F5_2)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * ewz/(sigx1 * wzdeno)
  }
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * prC/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1 + gv2, apply(F6_1, 1, sum) * ewz/(sigx1 * wzdeno),
    apply(F6_2, 1, sum) * prC/sigx1, sweep(Zvar, MARGIN = 1, STATS = prC * ewz *
      (sDiv1 - sDiv2)/(sigx1 * wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cauchit specification class membership
cgradmisfweibullnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewvsr <- exp(Wv/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  F1_1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k1)), MARGIN = 1, STATS = ewusr1,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F1_2 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k2)), MARGIN = 1, STATS = ewusr2,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"), 1, sum)
  sigx1 <- (ewz2 * sDiv2 + ewz1 * sDiv1)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1, STATS = ewusr1/ewvsr^3,
    FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1, STATS = ewusr2/ewvsr^3,
    FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1_1^2), MARGIN = 1, STATS = 1/ewvsr^2, FUN = "*") -
    0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1_2^2), MARGIN = 1, STATS = 1/ewvsr^2, FUN = "*") -
    0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 - FiMat)) * F1_1,
    MARGIN = 1, STATS = S * ewusr1/(k1^2 * ewvsr^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 - FiMat)) * F1_2,
    MARGIN = 1, STATS = S * ewusr2/(k2^2 * ewvsr^3), FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * ewz1/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * ewz2/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * ewz1/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * ewz2/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1 + gv2, apply(F6_1, 1, sum) * ewz1/sigx1, apply(F6_2,
    1, sum) * ewz2/sigx1, sweep(Zvar, MARGIN = 1, STATS = (sDiv1 - sDiv2)/(pi *
    sigx1 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# probit specification class membership
cgradmisfweibullnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewvsr <- exp(Wv/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  F1_1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k1)), MARGIN = 1, STATS = ewusr1,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F1_2 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k2)), MARGIN = 1, STATS = ewusr2,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"), 1, sum)
  sigx1 <- ((1 - pwZ) * sDiv2 + pwZ * sDiv1)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1, STATS = ewusr1/ewvsr^3,
    FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1, STATS = ewusr2/ewvsr^3,
    FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1_1^2), MARGIN = 1, STATS = 1/ewvsr^2, FUN = "*") -
    0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1_2^2), MARGIN = 1, STATS = 1/ewvsr^2, FUN = "*") -
    0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 - FiMat)) * F1_1,
    MARGIN = 1, STATS = S * ewusr1/(k1^2 * ewvsr^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 - FiMat)) * F1_2,
    MARGIN = 1, STATS = S * ewusr2/(k2^2 * ewvsr^3), FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pwZ/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * (1 - pwZ)/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pwZ/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * (1 - pwZ)/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1 + gv2, apply(F6_1, 1, sum) * pwZ/sigx1, apply(F6_2,
    1, sum) * (1 - pwZ)/sigx1, sweep(Zvar, MARGIN = 1, STATS = dwZ * (sDiv1 -
    sDiv2)/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cloglog specification class membership
cgradmisfweibullnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  k1 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 1]
  k2 <- parm[nXvar + 2 * nuZUvar + nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 3):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  F1_1 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k1)), MARGIN = 1, STATS = ewusr1,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F1_2 <- sweep(sweep((S * (-log(1 - FiMat))^(1/k2)), MARGIN = 1, STATS = ewusr2,
    FUN = "*"), MARGIN = 1, STATS = epsilon, FUN = "+")
  F2_1 <- dnorm(sweep(F1_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F2_2 <- dnorm(sweep(F1_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3_1 <- sweep(F2_1 * F1_1, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  F3_2 <- sweep(F2_2 * F1_2, MARGIN = 1, STATS = 1/ewvsr^3, FUN = "*")
  sDiv1 <- apply(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"), 1, sum)
  sDiv2 <- apply(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"), 1, sum)
  sigx1 <- ((1 - prZ) * sDiv1 + prZ * sDiv2)
  F4_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * F1_1, MARGIN = 1, STATS = ewusr1/ewvsr^3,
    FUN = "*")
  F4_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * F1_2, MARGIN = 1, STATS = ewusr2/ewvsr^3,
    FUN = "*")
  F5_1 <- sweep(sweep(0.5 * (F2_1 * F1_1^2), MARGIN = 1, STATS = 1/ewvsr^2, FUN = "*") -
    0.5 * F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  F5_2 <- sweep(sweep(0.5 * (F2_2 * F1_2^2), MARGIN = 1, STATS = 1/ewvsr^2, FUN = "*") -
    0.5 * F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*")
  F6_1 <- sweep((-log(1 - FiMat))^(1/k1) * F2_1 * log(-log(1 - FiMat)) * F1_1,
    MARGIN = 1, STATS = S * ewusr1/(k1^2 * ewvsr^3), FUN = "*")
  F6_2 <- sweep((-log(1 - FiMat))^(1/k2) * F2_2 * log(-log(1 - FiMat)) * F1_2,
    MARGIN = 1, STATS = S * ewusr2/(k2^2 * ewvsr^3), FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F3_1, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * (1 - prZ)/sigx1
    gx2[, k] <- apply(sweep(F3_2, MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * prZ/sigx1
  }
  gx <- gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F4_1)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F4_2)), MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * (1 - prZ)/sigx1
    gv2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * prZ/sigx1
  }
  gradll <- cbind(gx, gu1, gu2, gv1 + gv2, apply(F6_1, 1, sum) * (1 - prZ)/sigx1,
    apply(F6_2, 1, sum) * prZ/sigx1, sweep(Zvar, MARGIN = 1, STATS = prZ * ewz *
      (sDiv1 - sDiv2)/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for misf weibull-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
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
# logit specification class membership
misfweibullnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg, tol,
  gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar,
    Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  if (randStart) {
    rd <- rnorm(length(startVal), sd = sdStart)
    rd[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)] <- abs(rd[(nXvar +
      2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)])
    startVal <- startVal + rd
  }
  startLoglik <- sum(cmisfweibullnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfweibullnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfweibullnormlike_logit,
    grad = cgradmisfweibullnormlike_logit, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfweibullnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
        Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
        prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), control = list(iter.max = itermax, trace = printInfo,
      eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfweibullnormlike_logit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradmisfweibullnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmisfweibullnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfweibullnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initWeibull = initWeibull))
}

# cauchit specification class membership
misfweibullnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg, tol,
  gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar,
    Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  if (randStart) {
    rd <- rnorm(length(startVal), sd = sdStart)
    rd[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)] <- abs(rd[(nXvar +
      2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)])
    startVal <- startVal + rd
  }
  startLoglik <- sum(cmisfweibullnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfweibullnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfweibullnormlike_cauchit,
    grad = cgradmisfweibullnormlike_cauchit, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfweibullnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
        Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(cgradmisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
        prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisfweibullnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), control = list(iter.max = itermax, trace = printInfo,
      eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfweibullnormlike_cauchit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmisfweibullnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmisfweibullnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfweibullnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initWeibull = initWeibull))
}

# probit specification class membership
misfweibullnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg, tol,
  gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar,
    Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  if (randStart) {
    rd <- rnorm(length(startVal), sd = sdStart)
    rd[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)] <- abs(rd[(nXvar +
      2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)])
    startVal <- startVal + rd
  }
  startLoglik <- sum(cmisfweibullnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfweibullnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfweibullnormlike_probit,
    grad = cgradmisfweibullnormlike_probit, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfweibullnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
        Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(cgradmisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
        prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisfweibullnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), control = list(iter.max = itermax, trace = printInfo,
      eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfweibullnormlike_probit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmisfweibullnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmisfweibullnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfweibullnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initWeibull = initWeibull))
}

# cloglog specification class membership
misfweibullnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg, tol,
  gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar,
    Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  if (randStart) {
    rd <- rnorm(length(startVal), sd = sdStart)
    rd[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)] <- abs(rd[(nXvar +
      2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + 2)])
    startVal <- startVal + rd
  }
  startLoglik <- sum(cmisfweibullnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfweibullnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfweibullnormlike_cloglog,
    grad = cgradmisfweibullnormlike_cloglog, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfweibullnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
        Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(cgradmisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
        prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisfweibullnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      Zvar = Zvar, nZHvar = nZHvar)), control = list(iter.max = itermax, trace = printInfo,
      eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfweibullnormlike_cloglog(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmisfweibullnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmisfweibullnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfweibullnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initWeibull = initWeibull))
}

# Conditional efficiencies estimation ----------
#' efficiencies for misf weibull-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# logit specification class membership
cmisfweibullnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
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
    density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur1)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
      k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur2)/exp(Wv[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
      k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon2
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
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur1)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv[i]/2), k = k1, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur2)/exp(Wv[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv[i]/2), k = k2, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmisfweibullnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
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
    density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur1)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
      k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur2)/exp(Wv[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
      k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon2
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
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur1)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv[i]/2), k = k1, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur2)/exp(Wv[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv[i]/2), k = k2, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmisfweibullnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
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
    density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur1)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
      k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur2)/exp(Wv[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
      k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon2
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
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur1)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv[i]/2), k = k1, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur2)/exp(Wv[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv[i]/2), k = k2, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
cmisfweibullnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
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
    density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur1)/exp(Wv[i]/2)))
    u_c1[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
      k = k1, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon1
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur2)/exp(Wv[i]/2)))
    u_c2[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
      k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
      tol = 1e-15)$integral/density_epsilon2
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
      ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
      density_epsilon1 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur1)/exp(Wv[i]/2)))
      teBC_c1[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu1[i]/2),
        sigmaV = exp(Wv[i]/2), k = k1, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
      density_epsilon2 <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
        ur2)/exp(Wv[i]/2)))
      teBC_c2[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0, upperLimit = Inf,
        maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2), sigmaV = exp(Wv[i]/2),
        k = k2, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu2[i]/2),
        sigmaV = exp(Wv[i]/2), k = k2, epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
    }
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
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
#' marginal impact on efficiencies for misf weibull-normal distribution
#' @param object object of class sfacross
#' @noRd
# logit specification class membership
cmisfmargweibullnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu1/2) *
    gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu2/2) *
    gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargweibullnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmisfmargweibullnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu1/2) *
    gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu2/2) *
    gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargweibullnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmisfmargweibullnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu1/2) *
    gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu2/2) *
    gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargweibullnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmisfmargweibullnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu1/2) *
    gamma(1 + 1/k1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2, nrow = 1), matrix(exp(Wu2/2) *
    gamma(1 + 1/k2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargweibullnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  k1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 1]
  k2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar +
    2)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur1 <- exp(Wu1[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k1)
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur1)/exp(Wv[i]/2)))
    ur2 <- exp(Wu2[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k2)
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S * ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (gamma(1 + 2/k1) - (gamma(1 + 1/k1))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (gamma(1 + 2/k2) - (gamma(1 + 1/k2))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
