################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Model                           #
# Two types: - Common noise component (sigma_v)                                #
#            - Different noise component (multimodal noise - mnsf)             #
# Convolution: lognormal - normal                                              #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf lognormal-normal distribution
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
# Same sigma_v
czisflognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  for (i in 1:N) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(FiMat[i, ]))
    Pi1[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] + S *
      ur)/exp(Wv[i]/2)))
  }
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# Different sigma_v
cmnsflognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- numeric(N)
  for (i in 1:N) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(FiMat[i, ]))
    Pi1[i] <- mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      S * ur)/exp(Wv1[i]/2)))
  }
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf lognormal-normal distribution
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
# Same sigma_v
cstzisflognorm <- function(olsObj, epsiRes, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  Zvar, nZHvar, N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + lognormal - normal distributions...\n")
  initLog <- maxLik(logLik = clognormlike, start = cstlognorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE], nmuZUvar = 1, muHvar = muHvar[, 1,
      drop = FALSE]), grad = cgradlognormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    nmuZUvar = 1, muHvar = as.matrix(muHvar[, 1]), N = N,
    FiMat = FiMat)
  Esti <- initLog$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], rep(0, nmuZUvar -
    1), Esti[nXvar + 2], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 3], if (nvZVvar > 1) rep(0, nvZVvar -
    1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("ZI_", colnames(Zvar)))
  names(initLog$estimate) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)[1]), paste0("Zu_", colnames(uHvar)[1]),
    paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initLog = initLog))
}

# Different sigma_v
cstmnsflognorm <- function(olsObj, epsiRes, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  Zvar, nZHvar, N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + lognormal - normal distributions...\n")
  initLog <- maxLik(logLik = clognormlike, start = cstlognorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE], nmuZUvar = 1, muHvar = muHvar[, 1,
      drop = FALSE]), grad = cgradlognormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    nmuZUvar = 1, muHvar = as.matrix(muHvar[, 1]), N = N,
    FiMat = FiMat)
  Esti <- initLog$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], rep(0, nmuZUvar -
    1), Esti[nXvar + 2], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 3], if (nvZVvar > 1) rep(0, nvZVvar -
    1), Esti[nXvar + 3], if (nvZVvar > 1) rep(0, nvZVvar -
    1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("ZI_",
    colnames(Zvar)))
  names(initLog$estimate) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)[1]), paste0("Zu_", colnames(uHvar)[1]),
    paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initLog = initLog))
}

# Gradient of the likelihood function ----------
#' gradient for zisf lognormal-normal distribution
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
# Same sigma_v
cgradzisflognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
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
  qFi <- qnorm(FiMat)
  F1 <- exp(sweep(sweep(qFi, MARGIN = 1, STATS = ewu_h, FUN = "*"),
    MARGIN = 1, STATS = mu, FUN = "+"))
  F2 <- sweep((S * F1), MARGIN = 1, STATS = epsilon, FUN = "+")
  F3 <- sweep(F2, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F4 <- sweep(F2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  dF4 <- dnorm(F4, 0, 1)
  F5 <- dF4 * F3
  sDiv <- apply(sweep(dF4, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"),
    1, sum)
  sigx1 <- dnorm(S * (epsilon)/ewv_h, 0, 1)
  sigx2 <- (prC * sigx1/ewv_h + ewz * sDiv/(Q * wzdeno))
  sigx3 <- (1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2)
  sigx4 <- (sigx3 * sDiv - prC * sigx1/(wzdeno * ewv_h))
  sigx5 <- (Q * sigx2 * wzdeno)
  sigx6 <- (0.5 * (S^2 * sigx1 * (epsilon)^2/ewv_h^2) - 0.5 *
    sigx1)
  F6 <- sweep(dF4 * F1 * qFi * F3, MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F7 <- sweep((dF4 * F2^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*")
  F8 <- sweep((0.5 * F7 - 0.5 * dF4), MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  XF5 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF5[, k] <- apply(sweep(F5, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum) * ewz/(Q * wzdeno)/sigx2
  }
  gx <- XF5 + sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * sigx1 *
    (epsilon)/ewv_h^3/sigx2, FUN = "*")
  MUF <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    MUF[, k] <- apply(sweep(-(S * dF4 * F1 * F3), MARGIN = 1,
      STATS = muHvar[, k], FUN = "*"), 1, sum)
  }
  gmu <- sweep(MUF, MARGIN = 1, STATS = ewz/sigx5, FUN = "*")
  UF6 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF6[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
  }
  gu <- sweep(UF6, MARGIN = 1, STATS = ewz/sigx5, FUN = "*")
  VF8 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF8[, k] <- apply(sweep(F8, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  gv <- sweep(VF8, MARGIN = 1, STATS = ewz/(Q * wzdeno)/sigx2,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = sigx6 *
    prC/ewv_h/sigx2, FUN = "*")
  gradll <- cbind(gx, gmu, gu, gv, sweep(Zvar, MARGIN = 1,
    STATS = sigx4 * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v
cgradmnsflognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
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
  qFi <- qnorm(FiMat)
  F1 <- exp(sweep(sweep(qFi, MARGIN = 1, STATS = ewu_h, FUN = "*"),
    MARGIN = 1, STATS = mu, FUN = "+"))
  F2 <- sweep((S * F1), MARGIN = 1, STATS = epsilon, FUN = "+")
  F3 <- sweep(F2, MARGIN = 1, STATS = 1/ewv1_h^3, FUN = "*")
  dF4 <- dnorm(sweep(F2, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"),
    0, 1)
  F5 <- dF4 * F3
  sDiv <- apply(sweep(dF4, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"),
    1, sum)
  sigx1 <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  sigx2 <- (prC * sigx1/ewv2_h + ewz * sDiv/(Q * wzdeno))
  sigx3 <- (1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2)
  sigx4 <- (sigx3 * sDiv - prC * sigx1/(wzdeno * ewv2_h))
  sigx5 <- (Q * sigx2 * wzdeno)
  sigx6 <- (0.5 * (S^2 * sigx1 * (epsilon)^2/ewv2_h^2) - 0.5 *
    sigx1)
  F6 <- sweep(dF4 * F1 * qFi * F3, MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F7 <- sweep((dF4 * F2^2), MARGIN = 1, STATS = 1/ewv1_h^2,
    FUN = "*")
  F8 <- sweep((0.5 * F7 - 0.5 * dF4), MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F5, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = ewz/(Q * wzdeno), FUN = "*")
  XF3 <- XF2 + sweep(Xvar, MARGIN = 1, STATS = S^2 * prC *
    sigx1 * (epsilon)/ewv2_h^3, FUN = "*")
  gx <- sweep(XF3, MARGIN = 1, STATS = 1/sigx2, FUN = "*")
  MUF1 <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    MUF1[, k] <- apply(sweep(-(S * dF4 * F1 * F3), MARGIN = 1,
      STATS = muHvar[, k], FUN = "*"), 1, sum)
  }
  gmu <- sweep(MUF1, MARGIN = 1, STATS = ewz/sigx5, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
  }
  gu <- sweep(UF1, MARGIN = 1, STATS = ewz/sigx5, FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F8, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  gv1 <- sweep(VF1, MARGIN = 1, STATS = ewz/sigx5, FUN = "*")
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx6 * prC/(sigx2 *
    ewv2_h), FUN = "*")
  gz <- sweep(Zvar, MARGIN = 1, STATS = sigx4 * ewz/sigx2,
    FUN = "*")
  gradll <- cbind(gx, gmu, gu, gv1, gv2, gz)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf lognormal-normal distribution
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
# Same sigma_v
chesszisflognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
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
  qFi <- qnorm(FiMat)
  F1 <- exp(sweep(sweep(qFi, MARGIN = 1, STATS = ewu_h, FUN = "*"),
    MARGIN = 1, STATS = mu, FUN = "+"))
  F2 <- sweep((S * F1), MARGIN = 1, STATS = epsilon, FUN = "+")
  F3 <- sweep(F2, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  F4 <- sweep(F2, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  dF4 <- dnorm(F4, 0, 1)
  F5 <- sweep(S * F1 * F2, MARGIN = 1, STATS = 1/ewv_h^2, FUN = "*")
  sDiv <- apply(sweep(dF4, MARGIN = 1, STATS = 1/ewv_h, FUN = "*"),
    1, sum)
  sigx1 <- dnorm(S * (epsilon)/ewv_h, 0, 1)
  sigx2 <- (prC * sigx1/ewv_h + ewz * sDiv/(Q * wzdeno))
  sigx3 <- (1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2)
  sigx4 <- (sigx3 * sDiv - prC * sigx1/(wzdeno * ewv_h))
  sigx5 <- (Q * sigx2 * wzdeno)
  sigx6 <- (0.5 * (S^2 * sigx1 * (epsilon)^2/ewv_h^2) - 0.5 *
    sigx1)
  sigx7 <- (1/(Q * wzdeno) - (sigx4/sigx5 + Q/(Q * wzdeno)^2) *
    ewz)
  F6 <- sweep(dF4 * F1 * qFi * F3, MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F7 <- sweep((dF4 * F2^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*")
  F8 <- sweep((0.5 * F7 - 0.5 * dF4), MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  F9 <- sweep((F2^2), MARGIN = 1, STATS = 1/ewv_h^2, FUN = "*")
  F10 <- sweep((F9 - 1) * dF4, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  F11 <- sweep((F9 - 1) * dF4 * F1, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  F12 <- sweep((F9 - 1) * dF4 * F1 * qFi, MARGIN = 1, STATS = ewu_h/ewv_h^3,
    FUN = "*")
  F13 <- sweep(((1 - F5) * F2 + S * F1) * dF4 * F1, MARGIN = 1,
    STATS = 1/ewv_h^3, FUN = "*")
  F14 <- sweep(((1 - F5) * F2 + S * F1) * dF4 * F1 * qFi, MARGIN = 1,
    STATS = ewu_h/ewv_h^3, FUN = "*")
  F15 <- sweep((0.5 - 0.5 * (F5)) * qFi, MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F16 <- sweep((S * F1 * qFi), MARGIN = 1, STATS = ewu_h, FUN = "*")
  F17 <- sweep(F1 * qFi, MARGIN = 1, STATS = ewu_h/ewv_h^3,
    FUN = "*")
  F18 <- sweep(((0.5 * (0.5 * (F9) - 1) - 0.25) * dF4 * F9 -
    0.5 * (0.5 * F7 - 0.5 * dF4)), MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  XdF <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XdF[, k] <- apply(sweep(dF4 * F3, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  MUF <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    MUF[, k] <- apply(sweep(-(S * dF4 * F1 * F3), MARGIN = 1,
      STATS = muHvar[, k], FUN = "*"), 1, sum)
  }
  UF6 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF6[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
  }
  VF8 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF8[, k] <- apply(sweep(F8, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  Xsig1 <- sweep(XdF, MARGIN = 1, STATS = ewz/(Q * wzdeno),
    FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = S^2 * prC *
    sigx1 * (epsilon)/ewv_h^3, FUN = "*")
  VF9 <- sweep(VF8, MARGIN = 1, STATS = ewz/(Q * wzdeno), FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx6 * prC/ewv_h, FUN = "*")
  HX1 <- list()
  HXMU1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HMU1 <- list()
  HMUU1 <- list()
  HMUV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      F10[, r] * ewz/(Q * wzdeno)/sigx2, FUN = "*"), Xvar)
    HXMU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (S * F11)[, r] * ewz/sigx5, FUN = "*"), muHvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F12))[, r] * ewz/sigx5, FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      ((0.5 * (F9 - 2) - 0.5) * dF4 * F3)[, r] * ewz/(Q *
      wzdeno)/sigx2, FUN = "*"), vHvar)
    HMU1[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
      (S * F13)[, r] * ewz/sigx5, FUN = "*"), muHvar)
    HMUU1[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F14))[, r] * ewz/sigx5, FUN = "*"), uHvar)
    HMUV1[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
      (S * (0.5 + 0.5 * (2 - F9)) * dF4 * F1 * F3)[, r] *
      ewz/sigx5, FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (-(0.5 * (S * ((F15 + 0.5) * F2 + 0.5 * F16) * dF4 *
        F17)))[, r] * ewz/sigx5, FUN = "*"), uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (S * (0.25 + 0.5 * (1 - 0.5 * (F9))) * F6)[, r] *
      ewz/sigx5, FUN = "*"), vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      F18[, r] * ewz/(Q * wzdeno)/sigx2, FUN = "*"), vHvar)
  }
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar +
    nZHvar, ncol = nXvar + nmuZUvar + nuZUvar + nvZVvar +
    nZHvar)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * prC * sigx1 * (S^2 *
      (epsilon)^2/ewv_h^2 - 1)/ewv_h^3/sigx2, FUN = "*"),
    Xvar) - crossprod(sweep(Xsig1, MARGIN = 1, STATS = wHvar/sigx2^2,
    FUN = "*"), Xsig1)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- Reduce("+",
    HXMU1) - crossprod(sweep(Xsig1, MARGIN = 1, STATS = wHvar *
    ewz * Q * wzdeno/sigx5^2, FUN = "*"), MUF)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- Reduce("+", HXU1) - crossprod(sweep(Xsig1,
    MARGIN = 1, STATS = wHvar * ewz * Q * wzdeno/sigx5^2,
    FUN = "*"), UF6)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+", HXV1) +
    crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S^2 *
      (0.5 * (S^2 * (epsilon)^2/ewv_h^2 - 2) - 0.5) * prC *
      sigx1 * (epsilon)/ewv_h^3/sigx2, FUN = "*"), vHvar) -
    crossprod(sweep(Xsig1, MARGIN = 1, STATS = wHvar/sigx2^2,
      FUN = "*"), VF9)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(XdF,
    MARGIN = 1, STATS = wHvar * sigx3 * ewz/sigx2, FUN = "*") -
    (sweep(Xsig1, MARGIN = 1, STATS = wHvar * sigx4/sigx2 *
      ewz/sigx2, FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = wHvar *
      S^2 * prC * sigx1 * (epsilon)/(wzdeno * ewv_h^3) *
      ewz/sigx2, FUN = "*")), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar +
    nmuZUvar)] <- Reduce("+", HMU1) - crossprod(sweep(MUF,
    MARGIN = 1, STATS = wHvar * ewz^2/sigx5^2, FUN = "*"),
    MUF)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- Reduce("+", HMUU1) -
    crossprod(sweep(MUF, MARGIN = 1, STATS = wHvar * ewz^2/sigx5^2,
      FUN = "*"), UF6)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HMUV1) - crossprod(sweep(MUF, MARGIN = 1, STATS = wHvar *
    ewz/(sigx2 * sigx5), FUN = "*"), VF9)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar + nZHvar)] <- crossprod(sweep(MUF, MARGIN = 1,
    STATS = wHvar * sigx7 * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- Reduce("+",
    HU1) - crossprod(sweep(UF6, MARGIN = 1, STATS = wHvar *
    ewz^2/sigx5^2, FUN = "*"), UF6)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
      nuZUvar + nvZVvar)] <- Reduce("+", HUV1) - crossprod(sweep(UF6,
    MARGIN = 1, STATS = wHvar * ewz/(sigx2 * sigx5), FUN = "*"),
    VF9)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
      nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(UF6,
    MARGIN = 1, STATS = wHvar * sigx7 * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+", HV1) +
    crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar * prC *
      (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv_h^2) -
        1) - 0.25) * sigx1 * (epsilon)^2/ewv_h^2 - 0.5 *
        sigx6)/ewv_h/sigx2, FUN = "*"), vHvar) - crossprod(sweep(VF9,
    MARGIN = 1, STATS = wHvar/sigx2^2, FUN = "*"), VF9)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(VF8,
    MARGIN = 1, STATS = wHvar * sigx3 * ewz/sigx2, FUN = "*") -
    (sweep(VF9, MARGIN = 1, STATS = wHvar * sigx4/sigx2 *
      ewz/sigx2, FUN = "*") + sweep(vHvar, MARGIN = 1,
      STATS = wHvar * (0.5 * (S^2 * sigx1 * (epsilon)^2/(wzdeno *
        ewv_h^3)) - 0.5 * (wzdeno * sigx1 * ewv_h/(wzdeno *
        ewv_h)^2)) * prC * ewz/sigx2, FUN = "*")), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + nZHvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv_h) + ewv_h/(wzdeno *
      ewv_h)^2) * sigx1 - (sigx4^2/sigx2 + Q * (2 - 2 *
      (Q^2 * wzdeno * ewz/(Q * wzdeno)^2)) * sDiv/(Q *
      wzdeno)^2)) * ewz + sigx3 * sDiv - prC * sigx1/(wzdeno *
      ewv_h)) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Different sigma_v
chessmnsflognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar,
  nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
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
  qFi <- qnorm(FiMat)
  F1 <- exp(sweep(sweep(qFi, MARGIN = 1, STATS = ewu_h, FUN = "*"),
    MARGIN = 1, STATS = mu, FUN = "+"))
  F2 <- sweep((S * F1), MARGIN = 1, STATS = epsilon, FUN = "+")
  F3 <- sweep(F2, MARGIN = 1, STATS = 1/ewv1_h^3, FUN = "*")
  dF4 <- dnorm(sweep(F2, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"),
    0, 1)
  F5 <- dF4 * F3
  sDiv <- apply(sweep(dF4, MARGIN = 1, STATS = 1/ewv1_h, FUN = "*"),
    1, sum)
  sigx1 <- dnorm(S * (epsilon)/ewv2_h, 0, 1)
  sigx2 <- (prC * sigx1/ewv2_h + ewz * sDiv/(Q * wzdeno))
  sigx3 <- (1/(Q * wzdeno) - Q * ewz/(Q * wzdeno)^2)
  sigx4 <- (sigx3 * sDiv - prC * sigx1/(wzdeno * ewv2_h))
  sigx5 <- (Q * sigx2 * wzdeno)
  sigx6 <- (0.5 * (S^2 * sigx1 * (epsilon)^2/ewv2_h^2) - 0.5 *
    sigx1)
  F6 <- sweep(dF4 * F1 * qFi * F3, MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F7 <- sweep((dF4 * F2^2), MARGIN = 1, STATS = 1/ewv1_h^2,
    FUN = "*")
  F8 <- sweep((0.5 * F7 - 0.5 * dF4), MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(F5, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = ewz/(Q * wzdeno), FUN = "*")
  XF3 <- XF2 + sweep(Xvar, MARGIN = 1, STATS = S^2 * prC *
    sigx1 * (epsilon)/ewv2_h^3, FUN = "*")
  gx <- sweep(XF3, MARGIN = 1, STATS = 1/sigx2, FUN = "*")
  MUF1 <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    MUF1[, k] <- apply(sweep(-(S * dF4 * F1 * F3), MARGIN = 1,
      STATS = muHvar[, k], FUN = "*"), 1, sum)
  }
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep(-(0.5 * (S * F6)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
  }
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep(F8, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  sigx7 <- (S^2 * (epsilon)^2/ewv2_h^2 - 1)
  sigx8 <- (S^2 * (epsilon)^2/ewv2_h^2 - 2)
  sigx9 <- (1/(Q * wzdeno) - (sigx4/sigx5 + Q/(Q * wzdeno)^2) *
    ewz)
  sigx10 <- (0.5 * (0.5 * (S^2 * (epsilon)^2/ewv2_h^2) - 1) -
    0.25)
  sigx11 <- (S^2 * sigx10 * sigx1 * (epsilon)^2/(sigx2 * ewv2_h^3) -
    (sigx6 * prC + 0.5 * (sigx2 * ewv2_h)) * sigx6/(sigx2 *
      ewv2_h)^2)
  sigx12 <- ((sigx4 * sigx6/sigx2 + 0.5 * (S^2 * sigx1 * (epsilon)^2/(wzdeno *
    ewv2_h^2)))/ewv2_h - 0.5 * (wzdeno * sigx1 * ewv2_h/(wzdeno *
    ewv2_h)^2))
  sigx13 <- (2 - 2 * (Q^2 * wzdeno * ewz/(Q * wzdeno)^2))
  sigx14 <- ((prC * (1/(wzdeno^2 * ewv2_h) + ewv2_h/(wzdeno *
    ewv2_h)^2) * sigx1 - (sigx4^2/sigx2 + Q * sigx13 * sDiv/(Q *
    wzdeno)^2)) * ewz + sigx3 * sDiv - prC * sigx1/(wzdeno *
    ewv2_h))
  F9 <- sweep(F2^2, MARGIN = 1, STATS = 1/ewv1_h^2, FUN = "*")
  F10 <- sweep((F9 - 1) * dF4, MARGIN = 1, STATS = 1/ewv1_h^3,
    FUN = "*")
  F11 <- sweep(F2^2, MARGIN = 1, STATS = 1/ewv1_h^2, FUN = "*")
  F12 <- sweep((F11 - 1) * dF4 * F1, MARGIN = 1, STATS = 1/ewv1_h^3,
    FUN = "*")
  F13 <- sweep((F11 - 1) * dF4 * F1 * qFi, MARGIN = 1, STATS = ewu_h/ewv1_h^3,
    FUN = "*")
  F14 <- sweep(S * F1 * F2, MARGIN = 1, STATS = 1/ewv1_h^2,
    FUN = "*")
  F15 <- sweep(((1 - F14) * F2 + S * F1) * dF4 * F1 * qFi,
    MARGIN = 1, STATS = ewu_h/ewv1_h^3, FUN = "*")
  F16 <- sweep(((1 - F14) * F2 + S * F1) * dF4 * F1, MARGIN = 1,
    STATS = 1/ewv1_h^3, FUN = "*")
  F17 <- sweep((0.5 - 0.5 * (F14)) * qFi, MARGIN = 1, STATS = ewu_h,
    FUN = "*")
  F18 <- sweep((S * F1 * qFi), MARGIN = 1, STATS = ewu_h, FUN = "*")
  F19 <- sweep(((F17 + 0.5) * F2 + 0.5 * F18) * dF4 * F1 *
    qFi, MARGIN = 1, STATS = ewu_h/ewv1_h^3, FUN = "*")
  F20 <- sweep(((0.5 * (0.5 * (F11) - 1) - 0.25) * dF4 * F11 -
    0.5 * (0.5 * F7 - 0.5 * dF4)), MARGIN = 1, STATS = 1/ewv1_h,
    FUN = "*")
  HX1 <- list()
  HXMU1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HMU1 <- list()
  HMUU1 <- list()
  HMUV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      F10[, r] * ewz/(Q * wzdeno)/sigx2, FUN = "*"), Xvar)
    HXMU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (S * F12)[, r] * ewz/sigx5, FUN = "*"), muHvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F13))[, r] * ewz/sigx5, FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      ((0.5 * (F11 - 2) - 0.5) * F5)[, r] * ewz/sigx5,
      FUN = "*"), vHvar)
    HMU1[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
      (S * F16)[, r] * ewz/sigx5, FUN = "*"), muHvar)
    HMUU1[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * F15))[, r] * ewz/sigx5, FUN = "*"), uHvar)
    HMUV1[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
      (S * (0.5 + 0.5 * (2 - F11)) * dF4 * F1 * F3)[, r] *
      ewz/sigx5, FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (-(0.5 * (S * F19)))[, r] * ewz/sigx5, FUN = "*"),
      uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (S * (0.25 + 0.5 * (1 - 0.5 * (F11))) * F6)[, r] *
      ewz/sigx5, FUN = "*"), vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      (F20)[, r] * ewz/sigx5, FUN = "*"), vHvar)
  }
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + 2 *
    nvZVvar + nZHvar, ncol = nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * prC * sigx1 * sigx7/ewv2_h^3/sigx2,
    FUN = "*"), Xvar) - crossprod(sweep(XF3, MARGIN = 1,
    STATS = wHvar/sigx2^2, FUN = "*"), XF3)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- Reduce("+",
    HXMU1) - crossprod(sweep(XF3, MARGIN = 1, STATS = wHvar *
    ewz * Q * wzdeno/sigx5^2, FUN = "*"), MUF1)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- Reduce("+", HXU1) - crossprod(sweep(XF3,
    MARGIN = 1, STATS = wHvar * ewz * Q * wzdeno/sigx5^2,
    FUN = "*"), UF1)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+", HXV1) -
    crossprod(sweep(XF3, MARGIN = 1, STATS = wHvar * ewz *
      Q * wzdeno/sigx5^2, FUN = "*"), VF1)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * prC * S^2 * (0.5 * sigx8 -
      0.5) * sigx1 * (epsilon)/(sigx2 * ewv2_h^3), FUN = "*") -
    sweep(XF3, MARGIN = 1, STATS = wHvar * prC * sigx6 *
      ewv2_h/(sigx2 * ewv2_h)^2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(XF1,
    MARGIN = 1, STATS = wHvar * sigx3 * ewz/sigx2, FUN = "*") -
    (sweep(gx, MARGIN = 1, STATS = wHvar * sigx4 * ewz/sigx2,
      FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = wHvar *
      S^2 * prC * sigx1 * (epsilon)/(wzdeno * ewv2_h^3) *
      ewz/sigx2, FUN = "*")), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar +
    nmuZUvar)] <- Reduce("+", HMU1) - crossprod(sweep(MUF1,
    MARGIN = 1, STATS = wHvar * ewz * ewz/sigx5^2, FUN = "*"),
    MUF1)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- Reduce("+", HMUU1) -
    crossprod(sweep(MUF1, MARGIN = 1, STATS = wHvar * ewz^2/sigx5^2,
      FUN = "*"), UF1)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HMUV1) - crossprod(sweep(MUF1, MARGIN = 1, STATS = wHvar *
    ewz^2/sigx5^2, FUN = "*"), VF1)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(MUF1, MARGIN = 1, STATS = -wHvar *
    (sigx6 * prC * ewv2_h * ewz/(Q * (sigx2 * ewv2_h)^2 *
      wzdeno)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(MUF1, MARGIN = 1,
    STATS = wHvar * sigx9 * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- Reduce("+",
    HU1) - crossprod(sweep(UF1, MARGIN = 1, STATS = wHvar *
    ewz^2/sigx5^2, FUN = "*"), UF1)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
      nuZUvar + nvZVvar)] <- Reduce("+", HUV1) - crossprod(sweep(UF1,
    MARGIN = 1, STATS = wHvar * ewz^2/sigx5^2, FUN = "*"),
    VF1)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
      nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(UF1, MARGIN = 1,
    STATS = -wHvar * (sigx6 * prC * ewv2_h * ewz/(Q * (sigx2 *
      ewv2_h)^2 * wzdeno)), FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
      nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(UF1,
    MARGIN = 1, STATS = wHvar * sigx9 * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+", HV1) -
    crossprod(sweep(VF1, MARGIN = 1, STATS = wHvar * ewz^2/sigx5^2,
      FUN = "*"), VF1)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(VF1,
    MARGIN = 1, STATS = -wHvar * (sigx6 * prC * ewv2_h *
      ewz/(Q * (sigx2 * ewv2_h)^2 * wzdeno)), FUN = "*"),
    vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 2 *
    nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(VF1, MARGIN = 1, STATS = wHvar *
    sigx9 * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    prC * sigx11, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar), (nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx12 * prC * ewz/sigx2), FUN = "*"),
    Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar), (nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx14 * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf lognormal-normal distribution
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
# Same sigma_v
zisflognormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, muHvar,
  nmuZUvar, Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisflognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initLog <- start_st$initLog
  startVal <- start_st$StartVal
  startLoglik <- sum(czisflognormlike(startVal, nXvar = nXvar,
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
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisflognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
      FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisflognormlike,
    grad = cgradzisflognormlike, hess = chesszisflognormlike,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hs = function(parm) as(-chesszisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hess = function(parm) -chesszisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisflognormlike(mleObj$par,
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
      mleObj$hessian <- chesszisflognormlike(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisflognormlike(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisflognormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisflognormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initLog = initLog))
}

# Different sigma_v
mnsflognormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, muHvar,
  nmuZUvar, Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsflognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initLog <- start_st$initLog
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsflognormlike(startVal, nXvar = nXvar,
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
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsflognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
      FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmnsflognormlike,
    grad = cgradmnsflognormlike, hess = chessmnsflognormlike,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cmnsflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hs = function(parm) as(-chessmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(cmnsflognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hess = function(parm) -chessmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmnsflognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsflognormlike(mleObj$par,
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
      mleObj$hessian <- chessmnsflognormlike(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsflognormlike(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsflognormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsflognormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initLog = initLog))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf lognormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v
czisflognormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon <- (mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2))))
    u_c1[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), mu = mu[i], epsilon = epsilon[i],
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
    for (i in 1:object$Nobs) {
      ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
        ]))
      density_epsilon <- (mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv[i]/2))))
      teBC_c1[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
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
cmnsflognormeff <- function(object, level) {
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
    object$nuZUvar + 2 * object$nvZVvar + 1):(object$nXvar +
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv1[i]/2))))
    u_c1[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
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
    for (i in 1:object$Nobs) {
      ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
        ]))
      density_epsilon <- (mean(1/exp(Wv1[i]/2) * dnorm((epsilon[i] +
        object$S * ur)/exp(Wv1[i]/2))))
      teBC_c1[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv1[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal_c1[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv1[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
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
#' marginal impact on efficiencies for zisf lognormal-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v
czisfmarglognorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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

czisfmarglognorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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

# Different sigma_v
cmnsfmarglognorm_Eu <- function(object) {
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
    object$nuZUvar + 2 * object$nvZVvar + 1):(object$nXvar +
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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

cmnsfmarglognorm_Vu <- function(object) {
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
    object$nuZUvar + 2 * object$nvZVvar + 1):(object$nXvar +
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
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
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
