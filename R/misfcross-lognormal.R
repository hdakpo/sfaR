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
# Convolution: lognormal - normal                                              #
#------------------------------------------------------------------------------#
# Log-likelihood ----------
#' log-likelihood for misf lognormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
cmisflognormlike_logit <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  for (i in 1:N) {
    ur1 <- exp(mu1[i] + exp(Wu1[i]/2) * qnorm(FiMat[i, ]))
    ur2 <- exp(mu2[i] + exp(Wu2[i]/2) * qnorm(FiMat[i, ]))
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur1)/exp(Wv[i]/2)))
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# cauchit specification class membership
cmisflognormlike_cauchit <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  for (i in 1:N) {
    ur1 <- exp(mu1[i] + exp(Wu1[i]/2) * qnorm(FiMat[i, ]))
    ur2 <- exp(mu2[i] + exp(Wu2[i]/2) * qnorm(FiMat[i, ]))
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur1)/exp(Wv[i]/2)))
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# probit specification class membership
cmisflognormlike_probit <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  for (i in 1:N) {
    ur1 <- exp(mu1[i] + exp(Wu1[i]/2) * qnorm(FiMat[i, ]))
    ur2 <- exp(mu2[i] + exp(Wu2[i]/2) * qnorm(FiMat[i, ]))
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur1)/exp(Wv[i]/2)))
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# cloglog specification class membership
cmisflognormlike_cloglog <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  Pi2 <- numeric(N)
  for (i in 1:N) {
    ur1 <- exp(mu1[i] + exp(Wu1[i]/2) * qnorm(FiMat[i, ]))
    ur2 <- exp(mu2[i] + exp(Wu2[i]/2) * qnorm(FiMat[i, ]))
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur1)/exp(Wv[i]/2)))
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur2)/exp(Wv[i]/2)))
  }
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for misf lognormal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
cstmisflognorm <- function(olsObj, epsiRes, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  Zvar, nZHvar, N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + lognormal - normal distributions...\n")
  initLog <- maxLik(logLik = clognormlike, start = cstlognorm(olsObj = olsCoef,
    epsiRes = epsiRes, S = S, uHvar = as.matrix(uHvar[, 1]),
    nuZUvar = 1, vHvar = as.matrix(vHvar[, 1]), nvZVvar = 1,
    nmuZUvar = 1, muHvar = as.matrix(muHvar[, 1])), grad = cgradlognormlike,
    method = "BFGS", control = list(iterlim = itermax, printLevel = printInfo,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, nmuZUvar = 1,
    muHvar = as.matrix(muHvar[, 1]), N = N, FiMat = FiMat,
    wHvar = wHvar)
  Esti <- initLog$estimate
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nmuZUvar >
    1) rep(0, nmuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nmuZUvar >
    1) rep(0, nmuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 2], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 3], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)), paste0("Zmu_", colnames(muHvar)),
    paste0("Zu_", colnames(uHvar)), paste0("Zu_", colnames(uHvar)),
    paste0("Zv_", colnames(vHvar)), paste0("MI_", colnames(Zvar)))
  names(initLog$estimate) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)[1]), paste0("Zu_", colnames(uHvar)[1]),
    paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initLog = initLog))
}

# Gradient of the likelihood function ----------
#' gradient for misf lognormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
cgradmisflognormlike_logit <- function(parm, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv_h <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  F1_1 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewu1_h,
    FUN = "*"), MARGIN = 1, STATS = mu1, FUN = "+"))
  F1_2 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewu2_h,
    FUN = "*"), MARGIN = 1, STATS = mu2, FUN = "+"))
  F2_1 <- sweep(S * F1_1, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F2_2 <- sweep(S * F1_2, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F3_1 <- dnorm(sweep(F2_1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  F3_2 <- dnorm(sweep(F2_2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"))
  sDiv1 <- apply(sweep(F3_1, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"),
    1, sum)
  sDiv2 <- apply(sweep(F3_2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"),
    1, sum)
  F5_1 <- sweep(F3_1 * F2_1, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  F5_2 <- sweep(F3_2 * F2_2, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  sigx1 <- (prC * sDiv2 + ewz * sDiv1/wzdeno)
  F6_1 <- sweep(F3_1 * F1_1 * F2_1, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  F6_2 <- sweep(F3_2 * F1_2 * F2_2, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  F7_1 <- sweep(F3_1 * F1_1 * qnorm(FiMat) * F2_1, MARGIN = 1,
    STATS = ewu1_h/ewv_h^3, FUN = "*")
  F7_2 <- sweep(F3_2 * F1_2 * qnorm(FiMat) * F2_2, MARGIN = 1,
    STATS = ewu2_h/ewv_h^3, FUN = "*")
  F8_1 <- sweep(sweep(F3_1 * F2_1^2, MARGIN = 1, STATS = (0.5 *
    (1/ewv_h^2)), FUN = "*") - 0.5 * F3_1, MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  F8_2 <- sweep(sweep(F3_2 * F2_2^2, MARGIN = 1, STATS = (0.5 *
    (1/ewv_h^2)), FUN = "*") - 0.5 * F3_2, MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
  }
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gx <- gx1 + gx2
  gmu1 <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    gmu1[, k] <- apply(sweep(-(S * F6_1), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
  }
  gmu2 <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    gmu2[, k] <- apply(sweep(-(S * F6_2), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F7_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz/(sigx1 *
      wzdeno)
  }
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu2[, k] <- apply(sweep(-(0.5 * (S * F7_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F8_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz/(sigx1 * wzdeno)
  }
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv2[, k] <- apply(sweep(F8_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gradll <- cbind(gx, gmu1, gmu2, gu1, gu2, gv1 + gv2, sweep(Zvar,
    MARGIN = 1, STATS = prC * ewz * (sDiv1 - sDiv2)/(sigx1 *
      wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cauchit specification class membership
cgradmisflognormlike_cauchit <- function(parm, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
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
  F1_1 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewusr1,
    FUN = "*"), MARGIN = 1, STATS = mu1, FUN = "+"))
  F1_2 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewusr2,
    FUN = "*"), MARGIN = 1, STATS = mu2, FUN = "+"))
  F2_1 <- sweep(S * F1_1, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F2_2 <- sweep(S * F1_2, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F3_1 <- dnorm(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3_2 <- dnorm(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  sDiv1 <- apply(sweep(F3_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"),
    1, sum)
  sDiv2 <- apply(sweep(F3_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"),
    1, sum)
  F5_1 <- sweep(F3_1 * F2_1, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F5_2 <- sweep(F3_2 * F2_2, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  sigx1 <- (ewz2 * sDiv2 + ewz1 * sDiv1)
  F6_1 <- sweep(F3_1 * F1_1 * F2_1, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F6_2 <- sweep(F3_2 * F1_2 * F2_2, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F7_1 <- sweep(F3_1 * F1_1 * qnorm(FiMat) * F2_1, MARGIN = 1,
    STATS = ewusr1/ewvsr^3, FUN = "*")
  F7_2 <- sweep(F3_2 * F1_2 * qnorm(FiMat) * F2_2, MARGIN = 1,
    STATS = ewusr2/ewvsr^3, FUN = "*")
  F8_1 <- sweep(sweep(F3_1 * F2_1^2, MARGIN = 1, STATS = (0.5 *
    (1/ewvsr^2)), FUN = "*") - 0.5 * F3_1, MARGIN = 1, STATS = 1/ewvsr,
    FUN = "*")
  F8_2 <- sweep(sweep(F3_2 * F2_2^2, MARGIN = 1, STATS = (0.5 *
    (1/ewvsr^2)), FUN = "*") - 0.5 * F3_2, MARGIN = 1, STATS = 1/ewvsr,
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gx2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gx <- gx1 + gx2
  gmu1 <- matrix(nrow = N, ncol = nmuZUvar)
  gmu2 <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    gmu1[, k] <- apply(sweep(-(S * F6_1), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gmu2[, k] <- apply(sweep(-(S * F6_2), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F7_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz1/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F7_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F8_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz1/sigx1
    gv2[, k] <- apply(sweep(F8_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * ewz2/sigx1
  }
  gradll <- cbind(gx, gmu1, gmu2, gu1, gu2, gv1 + gv2, sweep(Zvar,
    MARGIN = 1, STATS = (sDiv1 - sDiv2)/(pi * sigx1 * ((Wz)^2 +
      1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# probit specification class membership
cgradmisflognormlike_probit <- function(parm, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
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
  F1_1 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewusr1,
    FUN = "*"), MARGIN = 1, STATS = mu1, FUN = "+"))
  F1_2 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewusr2,
    FUN = "*"), MARGIN = 1, STATS = mu2, FUN = "+"))
  F2_1 <- sweep(S * F1_1, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F2_2 <- sweep(S * F1_2, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F3_1 <- dnorm(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3_2 <- dnorm(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  sDiv1 <- apply(sweep(F3_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"),
    1, sum)
  sDiv2 <- apply(sweep(F3_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"),
    1, sum)
  F5_1 <- sweep(F3_1 * F2_1, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F5_2 <- sweep(F3_2 * F2_2, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  sigx1 <- ((1 - pwZ) * sDiv2 + pwZ * sDiv1)
  F6_1 <- sweep(F3_1 * F1_1 * F2_1, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F6_2 <- sweep(F3_2 * F1_2 * F2_2, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F7_1 <- sweep(F3_1 * F1_1 * qnorm(FiMat) * F2_1, MARGIN = 1,
    STATS = ewusr1/ewvsr^3, FUN = "*")
  F7_2 <- sweep(F3_2 * F1_2 * qnorm(FiMat) * F2_2, MARGIN = 1,
    STATS = ewusr2/ewvsr^3, FUN = "*")
  F8_1 <- sweep(sweep(F3_1 * F2_1^2, MARGIN = 1, STATS = (0.5 *
    (1/ewvsr^2)), FUN = "*") - 0.5 * F3_1, MARGIN = 1, STATS = 1/ewvsr,
    FUN = "*")
  F8_2 <- sweep(sweep(F3_2 * F2_2^2, MARGIN = 1, STATS = (0.5 *
    (1/ewvsr^2)), FUN = "*") - 0.5 * F3_2, MARGIN = 1, STATS = 1/ewvsr,
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gx2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gx <- gx1 + gx2
  gmu1 <- matrix(nrow = N, ncol = nmuZUvar)
  gmu2 <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    gmu1[, k] <- apply(sweep(-(S * F6_1), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gmu2[, k] <- apply(sweep(-(S * F6_2), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F7_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * pwZ/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F7_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F8_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * pwZ/sigx1
    gv2[, k] <- apply(sweep(F8_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - pwZ)/sigx1
  }
  gradll <- cbind(gx, gmu1, gmu2, gu1, gu2, gv1 + gv2, sweep(Zvar,
    MARGIN = 1, STATS = dwZ * (sDiv1 - sDiv2)/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cloglog specification class membership
cgradmisflognormlike_cloglog <- function(parm, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar +
    nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
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
  F1_1 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewusr1,
    FUN = "*"), MARGIN = 1, STATS = mu1, FUN = "+"))
  F1_2 <- exp(sweep(sweep(qnorm(FiMat), MARGIN = 1, STATS = ewusr2,
    FUN = "*"), MARGIN = 1, STATS = mu2, FUN = "+"))
  F2_1 <- sweep(S * F1_1, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F2_2 <- sweep(S * F1_2, MARGIN = 1, STATS = (epsilon), FUN = "+")
  F3_1 <- dnorm(sweep(F2_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  F3_2 <- dnorm(sweep(F2_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"))
  sDiv1 <- apply(sweep(F3_1, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"),
    1, sum)
  sDiv2 <- apply(sweep(F3_2, MARGIN = 1, STATS = 1/ewvsr, FUN = "*"),
    1, sum)
  F5_1 <- sweep(F3_1 * F2_1, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F5_2 <- sweep(F3_2 * F2_2, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  sigx1 <- ((1 - prZ) * sDiv1 + prZ * sDiv2)
  F6_1 <- sweep(F3_1 * F1_1 * F2_1, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F6_2 <- sweep(F3_2 * F1_2 * F2_2, MARGIN = 1, STATS = 1/ewvsr^3,
    FUN = "*")
  F7_1 <- sweep(F3_1 * F1_1 * qnorm(FiMat) * F2_1, MARGIN = 1,
    STATS = ewusr1/ewvsr^3, FUN = "*")
  F7_2 <- sweep(F3_2 * F1_2 * qnorm(FiMat) * F2_2, MARGIN = 1,
    STATS = ewusr2/ewvsr^3, FUN = "*")
  F8_1 <- sweep(sweep(F3_1 * F2_1^2, MARGIN = 1, STATS = (0.5 *
    (1/ewvsr^2)), FUN = "*") - 0.5 * F3_1, MARGIN = 1, STATS = 1/ewvsr,
    FUN = "*")
  F8_2 <- sweep(sweep(F3_2 * F2_2^2, MARGIN = 1, STATS = (0.5 *
    (1/ewvsr^2)), FUN = "*") - 0.5 * F3_2, MARGIN = 1, STATS = 1/ewvsr,
    FUN = "*")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(F5_1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gx2[, k] <- apply(sweep(F5_2, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gx <- gx1 + gx2
  gmu1 <- matrix(nrow = N, ncol = nmuZUvar)
  gmu2 <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    gmu1[, k] <- apply(sweep(-(S * F6_1), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gmu2[, k] <- apply(sweep(-(S * F6_2), MARGIN = 1, STATS = muHvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(-(0.5 * (S * F7_1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gu2[, k] <- apply(sweep(-(0.5 * (S * F7_2)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep(F8_1, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * (1 - prZ)/sigx1
    gv2[, k] <- apply(sweep(F8_2, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum) * prZ/sigx1
  }
  gradll <- cbind(gx, gmu1, gmu2, gu1, gu2, gv1 + gv2, sweep(Zvar,
    MARGIN = 1, STATS = prZ * ewz * (sDiv1 - sDiv2)/((1 -
      prZ) * sDiv1 + prZ * sDiv2), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for misf lognormal-normal distribution
#' @param start starting value for optimization
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
# Same sigma_u

## logit specification class membership
misflognormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  muHvar, nmuZUvar, Zvar, nZHvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisflognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initLog <- start_st$initLog
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisflognormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisflognormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisflognormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisflognormlike_logit,
    grad = cgradmisflognormlike_logit, hess = chessmisflognormlike_logit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hs = function(parm) as(-chessmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmisflognormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hess = function(parm) -chessmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisflognormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisflognormlike_logit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chessmisflognormlike_logit(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisflognormlike_logit(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisflognormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisflognormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initLog = initLog))
}

## cauchit specification class membership
misflognormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  muHvar, nmuZUvar, Zvar, nZHvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisflognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initLog <- start_st$initLog
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisflognormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisflognormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisflognormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisflognormlike_cauchit,
    grad = cgradmisflognormlike_cauchit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmisflognormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmisflognormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisflognormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisflognormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisflognormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisflognormlike_cauchit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmisflognormnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmisflognormnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmisflognormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisflognormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initLog = initLog))
}

## probit specification class membership
misflognormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  muHvar, nmuZUvar, Zvar, nZHvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisflognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initLog <- start_st$initLog
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisflognormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisflognormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisflognormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisflognormlike_probit,
    grad = cgradmisflognormlike_probit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmisflognormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmisflognormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisflognormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisflognormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisflognormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisflognormlike_probit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmisflognormnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmisflognormnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmisflognormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisflognormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initLog = initLog))
}

## cloglog specification class membership
misflognormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  muHvar, nmuZUvar, Zvar, nZHvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisflognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initLog <- start_st$initLog
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisflognormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisflognormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisflognormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisflognormlike_cloglog,
    grad = cgradmisflognormlike_cloglog, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmisflognormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisflognormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hs = as(jacobian(function(parm) -colSums(cgradmisflognormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmisflognormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisflognormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisflognormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmisflognormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisflognormlike_cloglog(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmisflognormnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradmisflognormnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmisflognormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisflognormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initLog = initLog))
}


# Conditional efficiencies estimation ----------
#' efficiencies for misf lognormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_u

## logit specification class membership
cmisflognormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
    Pi2[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + object$S *
      ur)/exp(Wv[i]/2)))
  }
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2))))
    density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2))))
    u_c1[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    u_c2[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
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
    for (i in 1:object$Nobs) {
      ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
        ]))
      density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2))))
      density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2))))
      teBC_c1[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_c2[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
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
cmisflognormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2))))
    density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2))))
    u_c1[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    u_c2[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
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
    for (i in 1:object$Nobs) {
      ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
        ]))
      density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2))))
      density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2))))
      teBC_c1[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_c2[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
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
cmisflognormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2))))
    density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2))))
    u_c1[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    u_c2[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
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
    for (i in 1:object$Nobs) {
      ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
        ]))
      density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2))))
      density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2))))
      teBC_c1[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_c2[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
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
cmisflognormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2))))
    density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv2[i]/2))))
    u_c1[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
    u_c2[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
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
    for (i in 1:object$Nobs) {
      ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
        ]))
      density_epsilon1 <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2))))
      density_epsilon2 <- (mean(1/exp(Wv2[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv2[i]/2))))
      teBC_c1[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_c2[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv2[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon1
      teBC_reciprocal_c2[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv2[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon2
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
#' marginal impact on efficiencies for misf lognormal-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_u

## logit specification class membership
cmisfmarglognorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2 + Wu)/2, ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmarglognorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(2 * (exp(Wu) - 1) * exp(2 * mu + exp(Wu)), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * exp(2 * mu + exp(Wu) + Wu), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmisfmarglognorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2 + Wu)/2, ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmarglognorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(2 * (exp(Wu) - 1) * exp(2 * mu + exp(Wu)), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * exp(2 * mu + exp(Wu) + Wu), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmisfmarglognorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2 + Wu)/2, ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmarglognorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(2 * (exp(Wu) - 1) * exp(2 * mu + exp(Wu)), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * exp(2 * mu + exp(Wu) + Wu), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmisfmarglognorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2 + Wu)/2, ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmarglognorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    2 * object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar +
    object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(object$Nobs)
  Pi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(2 * (exp(Wu) - 1) * exp(2 * mu + exp(Wu)), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * exp(2 * mu + exp(Wu) + Wu), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}
