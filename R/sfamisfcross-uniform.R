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
# Convolution: uniform - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for misf uniform-normal distribution
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
# logit specification class membership
cmisfuninormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# cauchit specification class membership
cmisfuninormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# probit specification class membership
cmisfuninormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# cloglog specification class membership
cmisfuninormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for misf uniform-normal distribution
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
#' @param @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik
#' @param tol parameter tolerance
#' @noRd
cstmisfuninorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstuninorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initUni <- NULL
  } else {
    cat("Initialization: SFA + uniform - normal distributions...\n")
    initUni <- maxLik::maxLik(logLik = cuninormlike, start = cstuninorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgraduninormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initUni$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1),
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("MISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initUni = initUni))
}

# Gradient of the likelihood function ----------
#' gradient for misf uniform-normal distribution
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
# logit specification class membership
cgradmisfuninormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
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
  musig1 <- (sqrt(12) * ewu1_h + S * epsilon)/ewv_h
  musig2 <- (sqrt(12) * ewu2_h + S * epsilon)/ewv_h
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewv_h
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- (prC * (depsi - dmusig2)/(sqrt(12) * ewu2_h) + (depsi - dmusig1) * ewz/(sqrt(12) *
    (wzdeno * ewu1_h)))
  sigx4 <- (prC * (pmusig2 - pepsi)/(sqrt(12) * ewu2_h) + ewz * (pmusig1 - pepsi)/(sqrt(12) *
    (wzdeno * ewu1_h)))
  sigx2 <- (sigx4 * ewv_h)
  sigx3 <- (dmusig1/(2 * (wzdeno * ewv_h)) - sqrt(12)/2 * (wzdeno * ewu1_h * (pmusig1 -
    pepsi)/(sqrt(12) * (wzdeno * ewu1_h))^2))
  sigx5 <- (dmusig2/(2 * ewv_h) - sqrt(12)/2 * (ewu2_h * (pmusig2 - pepsi)/(sqrt(12) *
    ewu2_h)^2))
  sigx6 <- (0.5 * (S * depsi * epsilon) - 0.5 * ((sqrt(12) * ewu1_h + S * epsilon) *
    dmusig1))
  sigx7 <- (0.5 * (S * depsi * epsilon) - 0.5 * ((sqrt(12) * ewu2_h + S * epsilon) *
    dmusig2))
  sigx8 <- (sigx6 * ewz/(sqrt(12) * (wzdeno * ewu1_h)) + sigx7 * prC/(sqrt(12) *
    ewu2_h))
  sigx9 <- (1/(sqrt(12) * (wzdeno * ewu1_h)) - sqrt(12) * (ewu1_h * ewz/(sqrt(12) *
    (wzdeno * ewu1_h))^2))
  sigx10 <- (sigx9 * (pmusig1 - pepsi) - prC * (pmusig2 - pepsi)/(sqrt(12) * (wzdeno *
    ewu2_h)))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx3 * ewz/sigx4, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = prC * sigx5/sigx4, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx8/sigx2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx10 *
      ewz/sigx4, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cauchit specification class membership
cgradmisfuninormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewvsr <- exp(Wv/2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (sqrt(12) * ewusr1 + S * (epsilon))/ewvsr
  musig2 <- (sqrt(12) * ewusr2 + S * (epsilon))/ewvsr
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewvsr
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- (ewz2 * (depsi - dmusig2)/(sqrt(12) * ewusr2) + ewz1 * (depsi - dmusig1)/(sqrt(12) *
    ewusr1))
  sigx2 <- (ewz2 * (pmusig2 - pepsi)/(sqrt(12) * ewusr2) + ewz1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1))
  sigx3_1 <- (dmusig1/(2 * ewvsr) - sqrt(12)/2 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1)^2))
  sigx3_2 <- (dmusig2/(2 * ewvsr) - sqrt(12)/2 * (ewusr2 * (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2)^2))
  sigx4_1 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr1 + S * (epsilon)) *
    dmusig1))
  sigx4_2 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr2 + S * (epsilon)) *
    dmusig2))
  sigx5 <- (ewz2 * sigx4_2/(sqrt(12) * ewusr2) + sigx4_1 * ewz1/(sqrt(12) * ewusr1))
  sigx6 <- ((pmusig1 - pepsi)/(sqrt(12) * ewusr1) - (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2))
  sigx7 <- (pi * sigx2 * ((Wz)^2 + 1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/(sigx2 * ewvsr), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ewz1 * sigx3_1/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = ewz2 * sigx3_2/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx5/(sigx2 * ewvsr), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx6/sigx7,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# probit specification class membership
cgradmisfuninormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewvsr <- exp(Wv/2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (sqrt(12) * ewusr1 + S * (epsilon))/ewvsr
  musig2 <- (sqrt(12) * ewusr2 + S * (epsilon))/ewvsr
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewvsr
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- ((1 - pwZ) * (depsi - dmusig2)/(sqrt(12) * ewusr2) + (depsi - dmusig1) *
    pwZ/(sqrt(12) * ewusr1))
  sigx2 <- ((1 - pwZ) * (pmusig2 - pepsi)/(sqrt(12) * ewusr2) + (pmusig1 - pepsi) *
    pwZ/(sqrt(12) * ewusr1))
  sigx3_1 <- (dmusig1/(2 * ewvsr) - sqrt(12)/2 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1)^2))
  sigx3_2 <- (dmusig2/(2 * ewvsr) - sqrt(12)/2 * (ewusr2 * (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2)^2))
  sigx4_1 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr1 + S * (epsilon)) *
    dmusig1))
  sigx4_2 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr2 + S * (epsilon)) *
    dmusig2))
  sigx5 <- (sigx4_1 * pwZ/(sqrt(12) * ewusr1) + sigx4_2 * (1 - pwZ)/(sqrt(12) *
    ewusr2))
  sigx6 <- ((pmusig1 - pepsi)/(sqrt(12) * ewusr1) - (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/(sigx2 * ewvsr), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = pwZ * sigx3_1/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = (1 - pwZ) * sigx3_2/sigx2, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx5/(sigx2 * ewvsr), FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx6 * dwZ/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cloglog specification class membership
cgradmisfuninormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewvsr <- exp(Wv/2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (sqrt(12) * ewusr1 + S * (epsilon))/ewvsr
  musig2 <- (sqrt(12) * ewusr2 + S * (epsilon))/ewvsr
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewvsr
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- ((1 - prZ) * (depsi - dmusig1)/(sqrt(12) * ewusr1) + (depsi - dmusig2) *
    prZ/(sqrt(12) * ewusr2))
  sigx2 <- ((1 - prZ) * (pmusig1 - pepsi)/(sqrt(12) * ewusr1) + prZ * (pmusig2 -
    pepsi)/(sqrt(12) * ewusr2))
  sigx3_1 <- (dmusig1/(2 * ewvsr) - sqrt(12)/2 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1)^2))
  sigx3_2 <- (dmusig2/(2 * ewvsr) - sqrt(12)/2 * (ewusr2 * (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2)^2))
  sigx4_1 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr1 + S * (epsilon)) *
    dmusig1))
  sigx4_2 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr2 + S * (epsilon)) *
    dmusig2))
  sigx5 <- (sigx4_1 * (1 - prZ)/(sqrt(12) * ewusr1) + sigx4_2 * prZ/(sqrt(12) *
    ewusr2))
  sigx6 <- ((pmusig1 - pepsi)/(sqrt(12) * ewusr1) - (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/(sigx2 * ewvsr), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (1 - prZ) * sigx3_1/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = prZ * sigx3_2/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx5/(sigx2 * ewvsr), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx6 *
      prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for misf uniform-normal distribution
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
# logit specification class membership
chessmisfuninormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
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
  musig1 <- (sqrt(12) * ewu1_h + S * epsilon)/ewv_h
  musig2 <- (sqrt(12) * ewu2_h + S * epsilon)/ewv_h
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewv_h
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- (prC * (depsi - dmusig2)/(sqrt(12) * ewu2_h) + (depsi - dmusig1) * ewz/(sqrt(12) *
    (wzdeno * ewu1_h)))
  sigx4 <- (prC * (pmusig2 - pepsi)/(sqrt(12) * ewu2_h) + ewz * (pmusig1 - pepsi)/(sqrt(12) *
    (wzdeno * ewu1_h)))
  sigx2 <- (sigx4 * ewv_h)
  sigx3 <- (dmusig1/(2 * (wzdeno * ewv_h)) - sqrt(12)/2 * (wzdeno * ewu1_h * (pmusig1 -
    pepsi)/(sqrt(12) * (wzdeno * ewu1_h))^2))
  sigx5 <- (dmusig2/(2 * ewv_h) - sqrt(12)/2 * (ewu2_h * (pmusig2 - pepsi)/(sqrt(12) *
    ewu2_h)^2))
  sigx6 <- (0.5 * (S * depsi * epsilon) - 0.5 * ((sqrt(12) * ewu1_h + S * epsilon) *
    dmusig1))
  sigx7 <- (0.5 * (S * depsi * epsilon) - 0.5 * ((sqrt(12) * ewu2_h + S * epsilon) *
    dmusig2))
  sigx8 <- (sigx6 * ewz/(sqrt(12) * (wzdeno * ewu1_h)) + sigx7 * prC/(sqrt(12) *
    ewu2_h))
  sigx9 <- (1/(sqrt(12) * (wzdeno * ewu1_h)) - sqrt(12) * (ewu1_h * ewz/(sqrt(12) *
    (wzdeno * ewu1_h))^2))
  sigx10 <- (sigx9 * (pmusig1 - pepsi) - prC * (pmusig2 - pepsi)/(sqrt(12) * (wzdeno *
    ewu2_h)))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((prC * (S * depsi * epsilon - (sqrt(12) * ewu2_h + S * epsilon) * dmusig2)/(sqrt(12) *
      ewu2_h) + ewz * (S * depsi * epsilon - (sqrt(12) * ewu1_h + S * epsilon) *
      dmusig1)/(sqrt(12) * (wzdeno * ewu1_h)))/(sigx4 * ewv_h^3) - sigx1^2/sigx2^2),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((sqrt(12) * ewu1_h + S * epsilon) * dmusig1/(2 * (wzdeno *
      ewv_h^2)) - (sigx1 * sigx3/sigx4 + sqrt(12)/2 * (wzdeno * (depsi - dmusig1) *
      ewu1_h/(sqrt(12) * (wzdeno * ewu1_h))^2))) * ewz/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sqrt(12) * ewu2_h + S * epsilon) * dmusig2/(2 *
      ewv_h^2) - (sigx1 * sigx5/sigx4 + sqrt(12)/2 * ((depsi - dmusig2) * ewu2_h/(sqrt(12) *
      ewu2_h)^2))) * prC/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((0.5 * (depsi * (S^2 * epsilon^2/ewv_h^2 -
      1)) - 0.5 * (((sqrt(12) * ewu1_h + S * epsilon)^2/ewv_h^2 - 1) * dmusig1)) *
      ewz/(sqrt(12) * (wzdeno * ewu1_h)) + (0.5 * (depsi * (S^2 * epsilon^2/ewv_h^2 -
      1)) - 0.5 * (((sqrt(12) * ewu2_h + S * epsilon)^2/ewv_h^2 - 1) * dmusig2)) *
      prC/(sqrt(12) * ewu2_h))/sigx2 - sigx8 * sigx1/sigx2^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx9 *
    (depsi - dmusig1) - (sigx1 * sigx10/sigx4 + prC * (depsi - dmusig2)/(sqrt(12) *
    (wzdeno * ewu2_h)))) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (((sqrt(12)/2 * (((sqrt(12)/2 * (dmusig1/ewv_h) -
      12 * (wzdeno^2 * ewu1_h * (pmusig1 - pepsi)/(sqrt(12) * (wzdeno * ewu1_h))^2)) *
      ewu1_h + 0.5 * (pmusig1 - pepsi)) * wzdeno/(sqrt(12) * (wzdeno * ewu1_h))^2) +
      sqrt(12)/2 * ((sqrt(12) * ewu1_h + S * epsilon) * dmusig1)/(2 * (wzdeno *
        ewv_h^3))) * ewu1_h + sigx3^2 * ewz/sigx4) * ewz/sigx4), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prC * sigx3 * sigx5 * ewz/sigx4^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx8 * sigx3 * ewv_h/sigx2^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 * ((sqrt(12) *
      ewu1_h + S * epsilon)^2/ewv_h^2)) * dmusig1)/(sqrt(12) * wzdeno) + sqrt(12)/2 *
      (sigx6 * wzdeno * ewu1_h/(sqrt(12) * (wzdeno * ewu1_h))^2))/sigx2) *
      ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((sqrt(12)/2 * (sigx9 * dmusig1/ewv_h) - (sqrt(12)/2 * wzdeno + sqrt(12) *
      ((0.5 - 12 * (wzdeno^2 * ewu1_h^2/(sqrt(12) * (wzdeno * ewu1_h))^2)) *
        ewz)) * (pmusig1 - pepsi)/(sqrt(12) * (wzdeno * ewu1_h))^2) * ewu1_h -
      sigx10 * sigx3 * ewz/sigx4) * ewz/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * ((prC *
    sigx5^2/sigx4 + (sqrt(12)/2 * (((sqrt(12)/2 * (dmusig2/ewv_h) - 12 * (ewu2_h *
    (pmusig2 - pepsi)/(sqrt(12) * ewu2_h)^2)) * ewu2_h + 0.5 * (pmusig2 - pepsi))/(sqrt(12) *
    ewu2_h)^2) + sqrt(12)/2 * ((sqrt(12) * ewu2_h + S * epsilon) * dmusig2)/(2 *
    ewv_h^3)) * ewu2_h) * prC/sigx4), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx8 * sigx5 * ewv_h/sigx2^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 * ((sqrt(12) *
      ewu2_h + S * epsilon)^2/ewv_h^2)) * dmusig2)/sqrt(12) + sqrt(12)/2 *
      (sigx7 * ewu2_h/(sqrt(12) * ewu2_h)^2))/sigx2) * prC), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((sigx10 * sigx5/sigx4 + dmusig2/(2 * (wzdeno * ewv_h)) -
      sqrt(12)/2 * (wzdeno * ewu2_h * (pmusig2 - pepsi)/(sqrt(12) * (wzdeno *
        ewu2_h))^2)) * prC * ewz/sigx4), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((0.25 * (S^3 * depsi * epsilon^3) - 0.25 * ((sqrt(12) *
      ewu1_h + S * epsilon)^3 * dmusig1)) * ewz/(sqrt(12) * (wzdeno * ewu1_h)) +
      (0.25 * (S^3 * depsi * epsilon^3) - 0.25 * ((sqrt(12) * ewu2_h + S *
        epsilon)^3 * dmusig2)) * prC/(sqrt(12) * ewu2_h))/(sigx4 * ewv_h^3) -
      sigx8 * (sigx6 * ewz/(sqrt(12) * (wzdeno * ewu1_h)) + sigx7 * prC/(sqrt(12) *
        ewu2_h) + 0.5 * sigx2)/sigx2^2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx6 * sigx9 - (sigx8 * sigx10/sigx4 + sigx7 *
      prC/(sqrt(12) * (wzdeno * ewu2_h)))) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(sqrt(12) * (wzdeno^2 * ewu2_h)) +
      sqrt(12) * (ewu2_h/(sqrt(12) * (wzdeno * ewu2_h))^2)) * (pmusig2 - pepsi) -
      (sigx10^2/sigx4 + (sqrt(12) * (1 - 24 * (wzdeno * ewu1_h^2 * ewz/(sqrt(12) *
        (wzdeno * ewu1_h))^2)) + sqrt(12)) * ewu1_h * (pmusig1 - pepsi)/(sqrt(12) *
        (wzdeno * ewu1_h))^2)) * ewz + sigx9 * (pmusig1 - pepsi) - prC *
      (pmusig2 - pepsi)/(sqrt(12) * (wzdeno * ewu2_h))) * ewz/sigx4, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cauchit specification class membership
chessmisfuninormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewvsr <- exp(Wv/2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (sqrt(12) * ewusr1 + S * (epsilon))/ewvsr
  musig2 <- (sqrt(12) * ewusr2 + S * (epsilon))/ewvsr
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewvsr
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- (ewz2 * (depsi - dmusig2)/(sqrt(12) * ewusr2) + ewz1 * (depsi - dmusig1)/(sqrt(12) *
    ewusr1))
  sigx2 <- (ewz2 * (pmusig2 - pepsi)/(sqrt(12) * ewusr2) + ewz1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1))
  sigx3_1 <- (dmusig1/(2 * ewvsr) - sqrt(12)/2 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1)^2))
  sigx3_2 <- (dmusig2/(2 * ewvsr) - sqrt(12)/2 * (ewusr2 * (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2)^2))
  sigx4_1 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr1 + S * (epsilon)) *
    dmusig1))
  sigx4_2 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr2 + S * (epsilon)) *
    dmusig2))
  sigx5 <- (ewz2 * sigx4_2/(sqrt(12) * ewusr2) + sigx4_1 * ewz1/(sqrt(12) * ewusr1))
  sigx6 <- ((pmusig1 - pepsi)/(sqrt(12) * ewusr1) - (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2))
  sigx7 <- (pi * sigx2 * ((Wz)^2 + 1))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((ewz2 * (S * depsi * (epsilon) - (sqrt(12) * ewusr2 + S * (epsilon)) * dmusig2)/(sqrt(12) *
      ewusr2) + ewz1 * (S * depsi * (epsilon) - (sqrt(12) * ewusr1 + S * (epsilon)) *
      dmusig1)/(sqrt(12) * ewusr1))/(sigx2 * ewvsr^3) - sigx1^2/(sigx2 * ewvsr)^2),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((sqrt(12) * ewusr1 + S * (epsilon)) * dmusig1/(2 * ewvsr^2) -
      (sigx1 * sigx3_1/sigx2 + sqrt(12)/2 * ((depsi - dmusig1) * ewusr1/(sqrt(12) *
        ewusr1)^2))) * ewz1/(sigx2 * ewvsr), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sqrt(12) * ewusr2 + S * (epsilon)) * dmusig2/(2 *
      ewvsr^2) - (sigx1 * sigx3_2/sigx2 + sqrt(12)/2 * ((depsi - dmusig2) *
      ewusr2/(sqrt(12) * ewusr2)^2))) * ewz2/(sigx2 * ewvsr), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((ewz2 * (0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 -
      1)) - 0.5 * (((sqrt(12) * ewusr2 + S * (epsilon))^2/ewvsr^2 - 1) * dmusig2))/(sqrt(12) *
      ewusr2) + (0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 - 1)) - 0.5 * (((sqrt(12) *
      ewusr1 + S * (epsilon))^2/ewvsr^2 - 1) * dmusig1)) * ewz1/(sqrt(12) *
      ewusr1))/(sigx2 * ewvsr) - sigx5 * sigx1/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (((depsi -
    dmusig1)/(sqrt(12) * ewusr1) - (depsi - dmusig2)/(sqrt(12) * ewusr2))/sigx7 -
    pi * sigx1 * sigx6 * ((Wz)^2 + 1)/sigx7^2)/ewvsr, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((ewz1 * sigx3_1^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 *
      (dmusig1/ewvsr) - 12 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) * ewusr1)^2)) *
      ewusr1 + 0.5 * (pmusig1 - pepsi))/(sqrt(12) * ewusr1)^2) + sqrt(12)/2 *
      ((sqrt(12) * ewusr1 + S * (epsilon)) * dmusig1)/(2 * ewvsr^3)) * ewusr1) *
      ewz1/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * ewz1 * sigx3_1 * sigx3_2/sigx2^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx5 * sigx3_1 * ewvsr/(sigx2 * ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewusr1 + S * (epsilon))^2/ewvsr^2)) * dmusig1)/sqrt(12) +
      sqrt(12)/2 * (sigx4_1 * ewusr1/(sqrt(12) * ewusr1)^2))/(sigx2 * ewvsr)) *
      ewz1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1/sigx7 - pi * sigx6 * ((Wz)^2 + 1) * ewz1/sigx7^2) * sigx3_1, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * ((ewz2 *
    sigx3_2^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 * (dmusig2/ewvsr) - 12 * (ewusr2 *
    (pmusig2 - pepsi)/(sqrt(12) * ewusr2)^2)) * ewusr2 + 0.5 * (pmusig2 - pepsi))/(sqrt(12) *
    ewusr2)^2) + sqrt(12)/2 * ((sqrt(12) * ewusr2 + S * (epsilon)) * dmusig2)/(2 *
    ewvsr^3)) * ewusr2) * ewz2/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx5 * sigx3_2 * ewvsr/(sigx2 * ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewusr2 + S * (epsilon))^2/ewvsr^2)) * dmusig2)/sqrt(12) +
      sqrt(12)/2 * (sigx4_2 * ewusr2/(sqrt(12) * ewusr2)^2))/(sigx2 * ewvsr)) *
      ewz2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((1/sigx7 + pi * sigx6 * ((Wz)^2 + 1) * ewz2/sigx7^2) *
      sigx3_2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((0.25 * (S^3 * depsi * (epsilon)^3) - 0.25 * ((sqrt(12) *
      ewusr1 + S * (epsilon))^3 * dmusig1)) * ewz1/(sqrt(12) * ewusr1) + (0.25 *
      (S^3 * depsi * (epsilon)^3) - 0.25 * ((sqrt(12) * ewusr2 + S * (epsilon))^3 *
      dmusig2)) * ewz2/(sqrt(12) * ewusr2))/(sigx2 * ewvsr^3) - sigx5 * (ewz2 *
      sigx4_2/(sqrt(12) * ewusr2) + sigx4_1 * ewz1/(sqrt(12) * ewusr1) + 0.5 *
      (sigx2 * ewvsr))/(sigx2 * ewvsr)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx4_1/(sqrt(12) * ewusr1) - sigx4_2/(sqrt(12) *
      ewusr2))/sigx7 - pi * sigx5 * sigx6 * ((Wz)^2 + 1)/sigx7^2)/ewvsr, FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx6 * ((pmusig1 - pepsi)/(sqrt(12) * ewusr1) +
      2 * (pi * Wz * sigx2) - (pmusig2 - pepsi)/(sqrt(12) * ewusr2))/sigx7^2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# probit specification class membership
chessmisfuninormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewvsr <- exp(Wv/2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (sqrt(12) * ewusr1 + S * (epsilon))/ewvsr
  musig2 <- (sqrt(12) * ewusr2 + S * (epsilon))/ewvsr
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewvsr
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- ((1 - pwZ) * (depsi - dmusig2)/(sqrt(12) * ewusr2) + (depsi - dmusig1) *
    pwZ/(sqrt(12) * ewusr1))
  sigx2 <- ((1 - pwZ) * (pmusig2 - pepsi)/(sqrt(12) * ewusr2) + (pmusig1 - pepsi) *
    pwZ/(sqrt(12) * ewusr1))
  sigx3_1 <- (dmusig1/(2 * ewvsr) - sqrt(12)/2 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1)^2))
  sigx3_2 <- (dmusig2/(2 * ewvsr) - sqrt(12)/2 * (ewusr2 * (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2)^2))
  sigx4_1 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr1 + S * (epsilon)) *
    dmusig1))
  sigx4_2 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr2 + S * (epsilon)) *
    dmusig2))
  sigx5 <- (sigx4_1 * pwZ/(sqrt(12) * ewusr1) + sigx4_2 * (1 - pwZ)/(sqrt(12) *
    ewusr2))
  sigx6 <- ((pmusig1 - pepsi)/(sqrt(12) * ewusr1) - (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (((1 - pwZ) * (S * depsi * (epsilon) - (sqrt(12) * ewusr2 + S * (epsilon)) *
      dmusig2)/(sqrt(12) * ewusr2) + pwZ * (S * depsi * (epsilon) - (sqrt(12) *
      ewusr1 + S * (epsilon)) * dmusig1)/(sqrt(12) * ewusr1))/(sigx2 * ewvsr^3) -
      sigx1^2/(sigx2 * ewvsr)^2), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((sqrt(12) * ewusr1 + S * (epsilon)) * dmusig1/(2 * ewvsr^2) -
      (sigx1 * sigx3_1/sigx2 + sqrt(12)/2 * ((depsi - dmusig1) * ewusr1/(sqrt(12) *
        ewusr1)^2))) * pwZ/(sigx2 * ewvsr), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sqrt(12) * ewusr2 + S * (epsilon)) * dmusig2/(2 *
      ewvsr^2) - (sigx1 * sigx3_2/sigx2 + sqrt(12)/2 * ((depsi - dmusig2) *
      ewusr2/(sqrt(12) * ewusr2)^2))) * (1 - pwZ)/(sigx2 * ewvsr), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((1 - pwZ) * (0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 -
      1)) - 0.5 * (((sqrt(12) * ewusr2 + S * (epsilon))^2/ewvsr^2 - 1) * dmusig2))/(sqrt(12) *
      ewusr2) + (0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 - 1)) - 0.5 * (((sqrt(12) *
      ewusr1 + S * (epsilon))^2/ewvsr^2 - 1) * dmusig1)) * pwZ/(sqrt(12) *
      ewusr1))/(sigx2 * ewvsr) - sigx5 * sigx1/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((depsi -
    dmusig1)/(sqrt(12) * ewusr1) - (sigx1 * sigx6/sigx2 + (depsi - dmusig2)/(sqrt(12) *
    ewusr2))) * dwZ/(sigx2 * ewvsr), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((pwZ * sigx3_1^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 *
      (dmusig1/ewvsr) - 12 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) * ewusr1)^2)) *
      ewusr1 + 0.5 * (pmusig1 - pepsi))/(sqrt(12) * ewusr1)^2) + sqrt(12)/2 *
      ((sqrt(12) * ewusr1 + S * (epsilon)) * dmusig1)/(2 * ewvsr^3)) * ewusr1) *
      pwZ/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * pwZ * sigx3_1 * sigx3_2/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx5 * sigx3_1 * ewvsr/(sigx2 * ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewusr1 + S * (epsilon))^2/ewvsr^2)) * dmusig1)/sqrt(12) +
      sqrt(12)/2 * (sigx4_1 * ewusr1/(sqrt(12) * ewusr1)^2))/(sigx2 * ewvsr)) *
      pwZ), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx6 * pwZ/sigx2) * sigx3_1 * dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (((1 -
    pwZ) * sigx3_2^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 * (dmusig2/ewvsr) -
    12 * (ewusr2 * (pmusig2 - pepsi)/(sqrt(12) * ewusr2)^2)) * ewusr2 + 0.5 *
    (pmusig2 - pepsi))/(sqrt(12) * ewusr2)^2) + sqrt(12)/2 * ((sqrt(12) * ewusr2 +
    S * (epsilon)) * dmusig2)/(2 * ewvsr^3)) * ewusr2) * (1 - pwZ)/sigx2), FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx5 * sigx3_2 * ewvsr/(sigx2 * ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewusr2 + S * (epsilon))^2/ewvsr^2)) * dmusig2)/sqrt(12) +
      sqrt(12)/2 * (sigx4_2 * ewusr2/(sqrt(12) * ewusr2)^2))/(sigx2 * ewvsr)) *
      (1 - pwZ)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((sigx6 * (1 - pwZ)/sigx2 + 1) * sigx3_2 * dwZ/sigx2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((0.25 * (S^3 * depsi * (epsilon)^3) - 0.25 * ((sqrt(12) *
      ewusr1 + S * (epsilon))^3 * dmusig1)) * pwZ/(sqrt(12) * ewusr1) + (0.25 *
      (S^3 * depsi * (epsilon)^3) - 0.25 * ((sqrt(12) * ewusr2 + S * (epsilon))^3 *
      dmusig2)) * (1 - pwZ)/(sqrt(12) * ewusr2))/(sigx2 * ewvsr^3) - sigx5 *
      ((1 - pwZ) * sigx4_2/(sqrt(12) * ewusr2) + sigx4_1 * pwZ/(sqrt(12) *
        ewusr1) + 0.5 * (sigx2 * ewvsr))/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx4_1/(sqrt(12) * ewusr1) - (sigx5 * sigx6/sigx2 +
      sigx4_2/(sqrt(12) * ewusr2))) * dwZ/(sigx2 * ewvsr), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx6 * dwZ/sigx2 + Wz) * sigx6 * dwZ/sigx2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cloglog specification class membership
chessmisfuninormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewvsr <- exp(Wv/2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (sqrt(12) * ewusr1 + S * (epsilon))/ewvsr
  musig2 <- (sqrt(12) * ewusr2 + S * (epsilon))/ewvsr
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi <- S * (epsilon)/ewvsr
  depsi <- dnorm(epsi, 0, 1)
  pepsi <- pnorm(epsi)
  sigx1 <- ((1 - prZ) * (depsi - dmusig1)/(sqrt(12) * ewusr1) + (depsi - dmusig2) *
    prZ/(sqrt(12) * ewusr2))
  sigx2 <- ((1 - prZ) * (pmusig1 - pepsi)/(sqrt(12) * ewusr1) + prZ * (pmusig2 -
    pepsi)/(sqrt(12) * ewusr2))
  sigx3_1 <- (dmusig1/(2 * ewvsr) - sqrt(12)/2 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) *
    ewusr1)^2))
  sigx3_2 <- (dmusig2/(2 * ewvsr) - sqrt(12)/2 * (ewusr2 * (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2)^2))
  sigx4_1 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr1 + S * (epsilon)) *
    dmusig1))
  sigx4_2 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * ((sqrt(12) * ewusr2 + S * (epsilon)) *
    dmusig2))
  sigx5 <- (sigx4_1 * (1 - prZ)/(sqrt(12) * ewusr1) + sigx4_2 * prZ/(sqrt(12) *
    ewusr2))
  sigx6 <- ((pmusig1 - pepsi)/(sqrt(12) * ewusr1) - (pmusig2 - pepsi)/(sqrt(12) *
    ewusr2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((prZ * (S * depsi * (epsilon) - (sqrt(12) * ewusr2 + S * (epsilon)) * dmusig2)/(sqrt(12) *
      ewusr2) + (1 - prZ) * (S * depsi * (epsilon) - (sqrt(12) * ewusr1 + S *
      (epsilon)) * dmusig1)/(sqrt(12) * ewusr1))/(sigx2 * ewvsr^3) - sigx1^2/(sigx2 *
      ewvsr)^2), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((sqrt(12) * ewusr1 + S * (epsilon)) * dmusig1/(2 * ewvsr^2) -
      (sigx1 * sigx3_1/sigx2 + sqrt(12)/2 * ((depsi - dmusig1) * ewusr1/(sqrt(12) *
        ewusr1)^2))) * (1 - prZ)/(sigx2 * ewvsr), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sqrt(12) * ewusr2 + S * (epsilon)) * dmusig2/(2 *
      ewvsr^2) - (sigx1 * sigx3_2/sigx2 + sqrt(12)/2 * ((depsi - dmusig2) *
      ewusr2/(sqrt(12) * ewusr2)^2))) * prZ/(sigx2 * ewvsr), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((prZ * (0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 -
      1)) - 0.5 * (((sqrt(12) * ewusr2 + S * (epsilon))^2/ewvsr^2 - 1) * dmusig2))/(sqrt(12) *
      ewusr2) + (0.5 * (depsi * (S^2 * (epsilon)^2/ewvsr^2 - 1)) - 0.5 * (((sqrt(12) *
      ewusr1 + S * (epsilon))^2/ewvsr^2 - 1) * dmusig1)) * (1 - prZ)/(sqrt(12) *
      ewusr1))/(sigx2 * ewvsr) - sigx5 * sigx1/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((depsi -
    dmusig1)/(sqrt(12) * ewusr1) - (sigx1 * sigx6/sigx2 + (depsi - dmusig2)/(sqrt(12) *
    ewusr2))) * prZ * ewz/(sigx2 * ewvsr), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (((1 - prZ) * sigx3_1^2/sigx2 + (sqrt(12)/2 *
      (((sqrt(12)/2 * (dmusig1/ewvsr) - 12 * (ewusr1 * (pmusig1 - pepsi)/(sqrt(12) *
        ewusr1)^2)) * ewusr1 + 0.5 * (pmusig1 - pepsi))/(sqrt(12) * ewusr1)^2) +
      sqrt(12)/2 * ((sqrt(12) * ewusr1 + S * (epsilon)) * dmusig1)/(2 * ewvsr^3)) *
      ewusr1) * (1 - prZ)/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * (1 - prZ) * sigx3_1 * sigx3_2/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx5 * sigx3_1 * ewvsr/(sigx2 * ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewusr1 + S * (epsilon))^2/ewvsr^2)) * dmusig1)/sqrt(12) +
      sqrt(12)/2 * (sigx4_1 * ewusr1/(sqrt(12) * ewusr1)^2))/(sigx2 * ewvsr)) *
      (1 - prZ)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx6 * (1 - prZ)/sigx2) * sigx3_1 * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * ((prZ *
    sigx3_2^2/sigx2 + (sqrt(12)/2 * (((sqrt(12)/2 * (dmusig2/ewvsr) - 12 * (ewusr2 *
    (pmusig2 - pepsi)/(sqrt(12) * ewusr2)^2)) * ewusr2 + 0.5 * (pmusig2 - pepsi))/(sqrt(12) *
    ewusr2)^2) + sqrt(12)/2 * ((sqrt(12) * ewusr2 + S * (epsilon)) * dmusig2)/(2 *
    ewvsr^3)) * ewusr2) * prZ/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx5 * sigx3_2 * ewvsr/(sigx2 * ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewusr2 + S * (epsilon))^2/ewvsr^2)) * dmusig2)/sqrt(12) +
      sqrt(12)/2 * (sigx4_2 * ewusr2/(sqrt(12) * ewusr2)^2))/(sigx2 * ewvsr)) *
      prZ), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((sigx6 * prZ/sigx2 + 1) * sigx3_2 * prZ * ewz/sigx2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((0.25 * (S^3 * depsi * (epsilon)^3) - 0.25 * ((sqrt(12) *
      ewusr1 + S * (epsilon))^3 * dmusig1)) * (1 - prZ)/(sqrt(12) * ewusr1) +
      (0.25 * (S^3 * depsi * (epsilon)^3) - 0.25 * ((sqrt(12) * ewusr2 + S *
        (epsilon))^3 * dmusig2)) * prZ/(sqrt(12) * ewusr2))/(sigx2 * ewvsr^3) -
      sigx5 * (prZ * sigx4_2/(sqrt(12) * ewusr2) + sigx4_1 * (1 - prZ)/(sqrt(12) *
        ewusr1) + 0.5 * (sigx2 * ewvsr))/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx4_1/(sqrt(12) * ewusr1) - (sigx5 * sigx6/sigx2 +
      sigx4_2/(sqrt(12) * ewusr2))) * prZ * ewz/(sigx2 * ewvsr), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx6 * (1 - (sigx6 * prZ/sigx2 + 1) * ewz) *
      prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for misf uniform-normal distribution
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
# logit specification class membership
misfuninormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfuninormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfuninormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfuninormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfuninormlike_logit,
    grad = cgradmisfuninormlike_logit, hess = chessmisfuninormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfuninormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfuninormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfuninormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfuninormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfuninormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfuninormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfuninormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfuninormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfuninormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfuninormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfuninormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfuninormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initUni = initUni))
}

# cauchit specification class membership
misfuninormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfuninormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfuninormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfuninormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfuninormlike_cauchit,
    grad = cgradmisfuninormlike_cauchit, hess = chessmisfuninormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfuninormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfuninormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfuninormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfuninormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfuninormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfuninormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfuninormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfuninormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfuninormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfuninormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfuninormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initUni = initUni))
}

# probit specification class membership
misfuninormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfuninormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfuninormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfuninormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfuninormlike_probit,
    grad = cgradmisfuninormlike_probit, hess = chessmisfuninormlike_probit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfuninormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfuninormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfuninormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfuninormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfuninormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfuninormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfuninormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfuninormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfuninormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfuninormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfuninormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initUni = initUni))
}

# cloglog specification class membership
misfuninormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfuninormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfuninormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfuninormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfuninormlike_cloglog,
    grad = cgradmisfuninormlike_cloglog, hess = chessmisfuninormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfuninormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfuninormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfuninormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfuninormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfuninormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfuninormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfuninormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfuninormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfuninormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfuninormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfuninormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfuninormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initUni = initUni))
}

# Conditional efficiencies estimation ----------
#' efficiencies for misf uniform-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# logit specification class membership
cmisfuninormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c2 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  theta_c1 <- ifelse(Group_c == 1, sqrt(12) * exp(Wu1/2), NA)
  theta_c2 <- ifelse(Group_c == 2, sqrt(12) * exp(Wu2/2), NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu1/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu2/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu1/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu2/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1, teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1, teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1, NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2, NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1, NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, teJLMS1_c = teJLMS1_c, teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c,
      teBC2_c = teBC2_c, teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      teBC1_c1 = teBC1_c1, teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c, teBC1_c2 = teBC1_c2,
      teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2, teBC2_reciprocal_c2 = teBC2_reciprocal_c2,
      ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1,
      effBC1_c2 = effBC1_c2, effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1,
      ReffBC2_c1 = ReffBC2_c1, ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2,
      theta_c1 = theta_c1, theta_c2 = theta_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1,
      u2_c1 = u2_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u1_c2 = u1_c2,
      u2_c2 = u2_c, ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, theta_c1 = theta_c1, theta_c2 = theta_c2)
  }
  return(res)
}

# cauchit specification class membership
cmisfuninormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c2 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  theta_c1 <- ifelse(Group_c == 1, sqrt(12) * exp(Wu1/2), NA)
  theta_c2 <- ifelse(Group_c == 2, sqrt(12) * exp(Wu2/2), NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu1/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu2/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu1/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu2/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1, teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1, teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1, NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2, NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1, NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, teJLMS1_c = teJLMS1_c, teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c,
      teBC2_c = teBC2_c, teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      teBC1_c1 = teBC1_c1, teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c, teBC1_c2 = teBC1_c2,
      teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2, teBC2_reciprocal_c2 = teBC2_reciprocal_c2,
      ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1,
      effBC1_c2 = effBC1_c2, effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1,
      ReffBC2_c1 = ReffBC2_c1, ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2,
      theta_c1 = theta_c1, theta_c2 = theta_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1,
      u2_c1 = u2_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u1_c2 = u1_c2,
      u2_c2 = u2_c, ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, theta_c1 = theta_c1, theta_c2 = theta_c2)
  }
  return(res)
}

# probit specification class membership
cmisfuninormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c2 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  theta_c1 <- ifelse(Group_c == 1, sqrt(12) * exp(Wu1/2), NA)
  theta_c2 <- ifelse(Group_c == 2, sqrt(12) * exp(Wu2/2), NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu1/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu2/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu1/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu2/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1, teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1, teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1, NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2, NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1, NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, teJLMS1_c = teJLMS1_c, teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c,
      teBC2_c = teBC2_c, teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      teBC1_c1 = teBC1_c1, teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c, teBC1_c2 = teBC1_c2,
      teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2, teBC2_reciprocal_c2 = teBC2_reciprocal_c2,
      ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1,
      effBC1_c2 = effBC1_c2, effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1,
      ReffBC2_c1 = ReffBC2_c1, ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2,
      theta_c1 = theta_c1, theta_c2 = theta_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1,
      u2_c1 = u2_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u1_c2 = u1_c2,
      u2_c2 = u2_c, ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, theta_c1 = theta_c1, theta_c2 = theta_c2)
  }
  return(res)
}

# cloglog specification class membership
cmisfuninormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- -exp(Wv/2) * ((dnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) - object$S *
    epsilon
  u2_c2 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 - pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  theta_c1 <- ifelse(Group_c == 1, sqrt(12) * exp(Wu1/2), NA)
  theta_c2 <- ifelse(Group_c == 2, sqrt(12) * exp(Wu2/2), NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu1/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S * epsilon +
      sqrt(12) * exp(Wu2/2))/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S * epsilon)/exp(Wv/2)) -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c2 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S * epsilon/exp(Wv/2) +
      exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu1/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu1/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu2/2))/exp(Wv/2) - exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((sqrt(12) * exp(Wu2/2) + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)))/(1 - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1, teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1, teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1, NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2, NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1, NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, teJLMS1_c = teJLMS1_c, teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c,
      teBC2_c = teBC2_c, teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      teBC1_c1 = teBC1_c1, teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c, teBC1_c2 = teBC1_c2,
      teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2, teBC2_reciprocal_c2 = teBC2_reciprocal_c2,
      ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1,
      effBC1_c2 = effBC1_c2, effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1,
      ReffBC2_c1 = ReffBC2_c1, ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2,
      theta_c1 = theta_c1, theta_c2 = theta_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u1_c = u1_c,
      u2_c = u2_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u1_c1 = u1_c1,
      u2_c1 = u2_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2, u1_c2 = u1_c2,
      u2_c2 = u2_c, ineff1_c1 = ineff1_c1, ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2,
      ineff2_c2 = ineff2_c2, theta_c1 = theta_c1, theta_c2 = theta_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for misf uniform-normal distribution
#' @param object object of class sfacross
#' @noRd
# logit specification class membership
cmisfmarguninorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarguninorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cauchit specification class membership
cmisfmarguninorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarguninorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# probit specification class membership
cmisfmarguninorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarguninorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cloglog specification class membership
cmisfmarguninorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(sqrt(3)/2 *
    exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarguninorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + 1):(object$nXvar +
    2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar + object$nZHvar)]
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
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) + object$S *
    epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
