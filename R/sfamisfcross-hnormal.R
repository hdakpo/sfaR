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
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for misf halfnormal-normal distribution
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
cmisfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# cauchit specification class membership
cmisfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# probit specification class membership
cmisfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# cloglog specification class membership
cmisfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu1) + exp(Wv))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu2) + exp(Wv))) *
    pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for misf halfnormal-normal distribution
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
cstmisfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA halfnormal - normal distribution...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1),
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("MISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for misf halfnormal-normal distribution
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
cgradmisfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewv <- exp(Wv)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv
  sigma_sq2 <- ewu2 + ewv
  sigmastar1 <- sqrt(ewu1 * ewv/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu1 * (epsilon)/ssq1)
  musig2 <- (S * ewu2 * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsi1 <- S * (epsilon)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon)/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sxq <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 * pmusig1 * (epsilon))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 * pmusig2 * (epsilon))
  wzsq1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  wzxsq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (2 * (prC * sigx1_2/sxq) + 2 * (sigx1_1 * ewz/wzsq1))
  sigx3_1 <- (depsi1 * ewz * pmusig1/wzxsq1)
  sigx3_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx4_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx4_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  usq1 <- (0.5 * (prU1 * ewv/sigmastar1) + sigmastar1)
  usq2 <- (0.5 * (prU2 * ewv/sigmastar2) + sigmastar2)
  sigx5_1 <- (1/ssq1 - usq1 * ewu1/ssq1^2)
  sigx5_2 <- (1/ssq2 - usq2 * ewu2/ssq2^2)
  sigx6_1 <- (0.5 * sigx4_1 - sigx5_1 * dmusig1 * depsi1)
  sigx6_2 <- (0.5 * sigx4_2 - sigx5_2 * dmusig2 * depsi2)
  s7sq2 <- (depsi2 * pmusig2/(sigma_sq2))
  sigx11_1 <- (wzdeno * depsi1 * pmusig1/wzxsq1^2)
  sigx11_2 <- (prC * depsi2 * pmusig2/(wzdeno * sqrt(sigma_sq2)))
  sigx7_1 <- (S * sigx6_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  sigx7_2 <- (S * sigx6_2 * (epsilon) - 0.5 * s7sq2)
  sigx8_1 <- ((2 * sigx3_2 + 2 * sigx3_1) * sqrt(sigma_sq1))
  sigx8_2 <- ((2 * sigx3_2 + 2 * sigx3_1) * sqrt(sigma_sq2))
  prV1 <- (1 - ewv/(sigma_sq1))
  prV2 <- (1 - ewv/(sigma_sq2))
  sigx9_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  s9sq2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx9_2 <- (s9sq2 * dmusig2 * depsi2 * ewu2/ssq2^2 + 0.5 * sigx4_2)
  sigx10_1 <- (sigx9_1 * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 * sigx4_1)
  sigx10_2 <- (S * sigx9_2 * (epsilon) - 0.5 * s7sq2)
  s12sq1 <- (1/wzxsq1 - ewz * sqrt(sigma_sq1)/wzxsq1^2)
  sigx12 <- (2 * (s12sq1 * depsi1 * pmusig1) - 2 * sigx11_2)
  sigx18_1 <- (S * sigx10_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/(2 * sigx3_2 + 2 *
    sigx3_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (ewu1 * ewz *
    sigx7_1/sigx8_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (prC *
    ewu2 * sigx7_2/sigx8_2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
    (ewv * ewz * sigx18_1/sigx8_1) + 2 * (prC * ewv * sigx10_2/sigx8_2), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx12 * ewz/(2 * sigx3_2 + 2 * sigx3_1),
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cauchit specification class membership
cgradmisfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/psq2 + ewz1 * depsi1 * pmusig1/psq1)
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx18 <- (pi * sigx2 * ((Wz)^2 + 1))
  sigx19 <- (sigx16_1 * ewz1/sqrt(sigma_sq1) + sigx16_2 * ewz2/sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_1 * ewz1/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = sigx11_2 * ewz2/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx19 * ewv/sigx2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9/sigx18,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# probit specification class membership
cgradmisfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/psq2 + depsi1 * pmusig1 * pwZ/psq1)
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx19 <- (sigx16_1 * pwZ/sqrt(sigma_sq1) + sigx16_2 * (1 - pwZ)/sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_1 * pwZ/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = sigx11_2 * (1 - pwZ)/sigx2, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx19 * ewv/sigx2, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx9 * dwZ/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cloglog specification class membership
cgradmisfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/psq1 + depsi2 * prZ * pmusig2/psq2)
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx19 <- (sigx16_1 * (1 - prZ)/sqrt(sigma_sq1) + sigx16_2 * prZ/sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_1 * (1 - prZ)/sigx2, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = sigx11_2 * prZ/sigx2, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx19 * ewv/sigx2, FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx9 * prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for misf halfnormal-normal distribution
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
chessmisfhalfnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv <- exp(Wv)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv
  sigma_sq2 <- ewu2 + ewv
  sigmastar1 <- sqrt(ewu1 * ewv/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  muvu1 <- (mu1 * ewv - S * ewu1 * (epsilon))
  muvu2 <- (mu2 * ewv - S * ewu2 * (epsilon))
  musig_sq1 <- muvu1/ssq1^2
  musig_sq2 <- muvu2/ssq2^2
  musig1 <- muvu1/ssq1
  musig2 <- muvu2/ssq2
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmu1 <- pnorm(mu1/ewu1_h)
  pmu2 <- pnorm(mu2/ewu2_h)
  dmu1 <- dnorm(mu1/ewu1_h, 0, 1)
  dmu2 <- dnorm(mu2/ewu2_h, 0, 1)
  musi1 <- (mu1 + S * (epsilon))
  musi2 <- (mu2 + S * (epsilon))
  epsi1 <- musi1/sqrt(sigma_sq1)
  epsi2 <- musi2/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * musi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * musi2 * pmusig2)
  spmu1 <- (wzdeno * (sigma_sq1) * pmu1 * sqrt(sigma_sq1))
  spmu2 <- ((sigma_sq2) * pmu2 * sqrt(sigma_sq2))
  sigx2 <- (prC * sigx1_2/spmu2 + sigx1_1 * ewz/spmu1)
  pmusq1 <- (wzdeno * pmu1 * sqrt(sigma_sq1))
  pmusq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx3 <- (prC * depsi2 * pmusig2/pmusq2 + depsi1 * ewz * pmusig1/pmusq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv/sigmastar1 - depsi1 * musi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv/sigmastar2 - depsi2 * musi2 * pmusig2)
  wup1 <- (pmusq1^2 * ewu1_h)
  wup2 <- (ewu2_h * pmusq2^2)
  sigx5_1 <- (sigx4_1/spmu1 - wzdeno * depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/wup1)
  sigx5_2 <- (sigx4_2/spmu2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/wup2)
  sigx6_1 <- (depsi1 * musi1^2 * pmusig1/(sigma_sq1)^2)
  sigx6_2 <- (depsi2 * musi2^2 * pmusig2/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  pvs1 <- (0.5 * (prU1 * ewv/sigmastar1) + sigmastar1)
  pvs2 <- (0.5 * (prU2 * ewv/sigmastar2) + sigmastar2)
  sigx7_1 <- (pvs1 * musig_sq1 + S * (epsilon)/ssq1)
  sigx7_2 <- (pvs2 * musig_sq2 + S * (epsilon)/ssq2)
  sigx8_1 <- (0.5 * sigx6_1 - sigx7_1 * dmusig1 * depsi1)
  sigx8_2 <- (0.5 * sigx6_2 - sigx7_2 * dmusig2 * depsi2)
  sigx9_1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  sigx9_2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx10_1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewu1_h)
  sigx10_2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewu2_h)
  sigx11_1 <- (sigx8_1 * ewu1/pmusq1 - (0.5 * sigx9_1 - 0.5 * sigx10_1) * wzdeno *
    depsi1 * pmusig1/pmusq1^2)
  sigx11_2 <- (sigx8_2 * ewu2/pmusq2 - (0.5 * sigx9_2 - 0.5 * sigx10_2) * depsi2 *
    pmusig2/pmusq2^2)
  prV1 <- (1 - ewv/(sigma_sq1))
  prV2 <- (1 - ewv/(sigma_sq2))
  sigx12_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/ssq1 - sigx12_1 * musig_sq1)
  sigx13_2 <- (mu2/ssq2 - sigx12_2 * musig_sq2)
  sigx14_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx14_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx15_1 <- (wzdeno * depsi1 * pmusig1 * pmu1/pmusq1^2)
  sigx15_2 <- (depsi2 * pmusig2 * pmu2/pmusq2^2)
  sigx16_1 <- (sigx14_1/(wzdeno * pmu1) - 0.5 * sigx15_1)
  sigx16_2 <- (sigx14_2/pmu2 - 0.5 * sigx15_2)
  sigx17 <- (1/pmusq1 - ewz * pmu1 * sqrt(sigma_sq1)/pmusq1^2)
  sigx18 <- (wzdeno * pmu2 * sqrt(sigma_sq2))
  sigx19 <- (sigx17 * depsi1 * pmusig1 - prC * depsi2 * pmusig2/sigx18)
  sigx20_1 <- ((musi1^2/(sigma_sq1) - 1) * pmusig1 + dmusig1 * ewu1 * musi1/ssq1)
  sigx20_2 <- ((musi2^2/(sigma_sq2) - 1) * pmusig2 + dmusig2 * ewu2 * musi2/ssq2)
  sigx21_1 <- (depsi1 * musi1 - depsi1 * muvu1/ewv)
  sigx21_2 <- (depsi2 * musi2 - depsi2 * muvu2/ewv)
  musq1 <- (musi1^2/(sigma_sq1) - 2)
  musq2 <- (musi2^2/(sigma_sq2) - 2)
  sigx22_1 <- ((musq1 * pmusig1 + dmusig1 * ewu1 * musi1/ssq1) * depsi1 * musi1/(sigma_sq1)^2)
  sigx22_2 <- ((musq2 * pmusig2 + dmusig2 * ewu2 * musi2/ssq2) * depsi2 * musi2/(sigma_sq2)^2)
  sigx23_1 <- sigx7_1 * depsi1 * musi1/(sigma_sq1)
  sigx23_2 <- sigx7_2 * depsi2 * musi2/(sigma_sq2)
  sigx24_1 <- (depsi1 * muvu1/ewu1 + depsi1 * musi1)
  sigx24_2 <- (depsi2 * muvu2/ewu2 + depsi2 * musi2)
  sigx25_1 <- (wzdeno^2 * (sigma_sq1) * pmu1^2/pmusq1^2)
  sigx26_1 <- wzdeno^2 * pmu1^2 * sqrt(sigma_sq1)/pmusq1^2
  sigx27_1 <- (((2 - musi1^2/(sigma_sq1)) * pmusig1 + dmusig1 * ewv * musi1/ssq1) *
    depsi1 * musi1/(sigma_sq1)^2)
  sigx27_2 <- (((2 - musi2^2/(sigma_sq2)) * pmusig2 + dmusig2 * ewv * musi2/ssq2) *
    depsi2 * musi2/(sigma_sq2)^2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * ((sigx20_1 * depsi1 + dmusig1 * sigx21_1 * ewu1/ssq1) * ewz/spmu1 +
    (sigx20_2 * depsi2 + dmusig2 * sigx21_2 * ewu2/ssq2) * prC/spmu2 - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((0.5 * sigx22_1 - (sigx23_1 + (pvs1 * ewu1/ssq1^2 -
      (sigx7_1 * muvu1/ewv + 1/sigmastar1)/(sigma_sq1)) * depsi1) * dmusig1) *
      ewu1/pmusq1 - (sigx11_1 * sigx2/sigx3 + (0.5 * sigx9_1 - 0.5 * sigx10_1) *
      wzdeno * sigx1_1/(pmusq1^2 * (sigma_sq1)))) * ewz/sigx3, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * sigx22_2 - (sigx23_2 + (pvs2 * ewu2/ssq2^2 -
      (sigx7_2 * muvu2/ewv + 1/sigmastar2)/(sigma_sq2)) * depsi2) * dmusig2) *
      ewu2/pmusq2 - (sigx11_2 * sigx2/sigx3 + (0.5 * sigx9_2 - 0.5 * sigx10_2) *
      sigx1_2/((sigma_sq2) * pmusq2^2))) * prC/sigx3, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((sigx21_1 * sigx13_1/(sigma_sq1) - sigx12_1 *
      depsi1 * ewu1/ssq1^2) * dmusig1 + 0.5 * sigx22_1)/(wzdeno * pmu1) - 0.5 *
      (wzdeno * sigx1_1 * pmu1/(pmusq1^2 * (sigma_sq1)))) * ewz/sqrt(sigma_sq1) +
      (((sigx21_2 * sigx13_2/(sigma_sq2) - sigx12_2 * depsi2 * ewu2/ssq2^2) *
        dmusig2 + 0.5 * sigx22_2)/pmu2 - 0.5 * (sigx1_2 * pmu2/((sigma_sq2) *
        pmusq2^2))) * prC/sqrt(sigma_sq2) - (sigx16_1 * ewz/sqrt(sigma_sq1) +
      sigx16_2 * prC/sqrt(sigma_sq2)) * sigx2/sigx3) * ewv/sigx3, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx17 *
    sigx1_1/(sigma_sq1) - (sigx2 * sigx19/sigx3 + prC * sigx1_2/(wzdeno * (sigma_sq2) *
    pmu2 * sqrt(sigma_sq2)))) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (musi1^2/(sigma_sq1)) - 2) *
      pmusig1/(sigma_sq1) - sigx7_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      musi1^2/(sigma_sq1)^2 - ((sigx7_1^2 * ewu1 * musig1 + ((0.5 * (ewu1/(sigma_sq1)) -
      0.5 * (0.5 * prU1 + ewu1/(sigma_sq1))) * prU1 * ewv * muvu1/sigmastar1 -
      pvs1 * (2 * (pvs1 * (sigma_sq1) * muvu1 * sigmastar1/ssq1^2) + 2 * (S *
        (epsilon))) * ewu1)/ssq1^2) * depsi1 + sigx7_1 * (0.5 * (depsi1 *
      ewu1 * musi1^2/(sigma_sq1)^2) + depsi1)) * dmusig1)/pmusq1 - sigx8_1 *
      (0.5 * sigx9_1 - 0.5 * sigx10_1) * wzdeno/pmusq1^2) * ewu1 - ((((0.5 *
      (((1 - 0.5 * (ewu1/(sigma_sq1))) * pmu1 - 0.5 * (mu1 * dmu1/ewu1_h)) *
        ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 * ((0.5 * (mu1^2/ewu1_h^2) - 0.5) *
      sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) * dmu1/ewu1_h)) * depsi1 +
      0.5 * ((0.5 * sigx9_1 - 0.5 * sigx10_1) * depsi1 * ewu1 * musi1^2/(sigma_sq1)^2)) *
      pmusig1 - (sigx7_1 * dmusig1 * ewu1 + 2 * ((0.5 * sigx9_1 - 0.5 * sigx10_1) *
      wzdeno^2 * pmusig1 * pmu1 * sqrt(sigma_sq1)/pmusq1^2)) * (0.5 * sigx9_1 -
      0.5 * sigx10_1) * depsi1) * wzdeno/pmusq1^2 + sigx11_1^2 * ewz/sigx3)) *
      ewz/sigx3, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * prC * ewz/sigx3^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx7_1 * depsi1 * muvu1/sigmastar1 + 0.5 * (depsi1 * musi1^2/(sigma_sq1))) *
      sigx13_1/(sigma_sq1) - ((0.5 * (prU1 * ewv/(sigma_sq1)) + 0.5 * ((ewu1/(sigma_sq1) -
      1) * ewv/(sigma_sq1) + 1 - 0.5 * (prU1 * prV1))) * muvu1/sigmastar1 +
      mu1 * pvs1 - sigx12_1 * (2 * (pvs1 * (sigma_sq1) * muvu1 * sigmastar1/ssq1^2) +
      S * (epsilon))) * depsi1/ssq1^2) * dmusig1 + 0.5 * (((0.5 * (musi1^2/(sigma_sq1)) -
      2) * pmusig1/(sigma_sq1) - sigx7_1 * dmusig1) * depsi1 * musi1^2/(sigma_sq1)^2))/(wzdeno *
      pmu1) - 0.5 * (sigx16_1/(sigma_sq1))) * ewu1 + (0.5 * (mu1 * sigx14_1 *
      dmu1/((wzdeno * pmu1)^2 * ewu1_h)) - 0.5 * ((sigx8_1 * ewu1 * pmu1 -
      (0.5 * (mu1 * dmu1/ewu1_h) + 2 * ((0.5 * sigx9_1 - 0.5 * sigx10_1) *
        sigx26_1)) * depsi1 * pmusig1)/pmusq1^2)) * wzdeno)/sqrt(sigma_sq1) -
      (sigx16_1 * ewz/sqrt(sigma_sq1) + sigx16_2 * prC/sqrt(sigma_sq2)) * sigx11_1/sigx3) *
    ewv * ewz/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (sigx17 * depsi1 * ewu1 * musi1^2/(sigma_sq1)^2) - ((2 - 2 * sigx25_1) *
      ewz + 1) * (0.5 * sigx9_1 - 0.5 * sigx10_1) * depsi1/pmusq1^2) * pmusig1 -
      (sigx7_1 * sigx17 * dmusig1 * depsi1 * ewu1 + sigx11_1 * sigx19 * ewz/sigx3)) *
    ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((0.5 *
    (((0.5 * (musi2^2/(sigma_sq2)) - 2) * pmusig2/(sigma_sq2) - sigx7_2 * dmusig2) *
      ewu2) + 0.5 * pmusig2) * depsi2 * musi2^2/(sigma_sq2)^2 - ((sigx7_2^2 *
    ewu2 * musig2 + ((0.5 * (ewu2/(sigma_sq2)) - 0.5 * (0.5 * prU2 + ewu2/(sigma_sq2))) *
    prU2 * ewv * muvu2/sigmastar2 - pvs2 * (2 * (pvs2 * (sigma_sq2) * muvu2 *
    sigmastar2/ssq2^2) + 2 * (S * (epsilon))) * ewu2)/ssq2^2) * depsi2 + sigx7_2 *
    (0.5 * (depsi2 * ewu2 * musi2^2/(sigma_sq2)^2) + depsi2)) * dmusig2)/pmusq2 -
    sigx8_2 * (0.5 * sigx9_2 - 0.5 * sigx10_2)/pmusq2^2) * ewu2 - ((((0.5 * (((1 -
    0.5 * (ewu2/(sigma_sq2))) * pmu2 - 0.5 * (mu2 * dmu2/ewu2_h)) * ewu2/sqrt(sigma_sq2)) -
    0.5 * (mu2 * ((0.5 * (mu2^2/ewu2_h^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewu2_h)) * depsi2 + 0.5 * ((0.5 * sigx9_2 - 0.5 * sigx10_2) * depsi2 *
    ewu2 * musi2^2/(sigma_sq2)^2)) * pmusig2 - (sigx7_2 * dmusig2 * ewu2 + 2 *
    ((0.5 * sigx9_2 - 0.5 * sigx10_2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/pmusq2^2)) *
    (0.5 * sigx9_2 - 0.5 * sigx10_2) * depsi2)/pmusq2^2 + sigx11_2^2 * prC/sigx3)) *
    prC/sigx3, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx7_2 * depsi2 * muvu2/sigmastar2 + 0.5 * (depsi2 * musi2^2/(sigma_sq2))) *
      sigx13_2/(sigma_sq2) - ((0.5 * (prU2 * ewv/(sigma_sq2)) + 0.5 * ((ewu2/(sigma_sq2) -
      1) * ewv/(sigma_sq2) + 1 - 0.5 * (prU2 * prV2))) * muvu2/sigmastar2 +
      mu2 * pvs2 - sigx12_2 * (2 * (pvs2 * (sigma_sq2) * muvu2 * sigmastar2/ssq2^2) +
      S * (epsilon))) * depsi2/ssq2^2) * dmusig2 + 0.5 * (((0.5 * (musi2^2/(sigma_sq2)) -
      2) * pmusig2/(sigma_sq2) - sigx7_2 * dmusig2) * depsi2 * musi2^2/(sigma_sq2)^2)) *
      ewu2 + 0.5 * (mu2 * sigx14_2 * dmu2/(ewu2_h * pmu2)))/pmu2 - (0.5 * ((sigx8_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewu2_h) + 2 * ((0.5 * sigx9_2 - 0.5 *
      sigx10_2) * pmu2^2 * sqrt(sigma_sq2)/pmusq2^2)) * depsi2 * pmusig2)/pmusq2^2) +
      0.5 * (sigx16_2 * ewu2/(sigma_sq2))))/sqrt(sigma_sq2) - (sigx16_1 * ewz/sqrt(sigma_sq1) +
      sigx16_2 * prC/sqrt(sigma_sq2)) * sigx11_2/sigx3) * prC * ewv/sigx3,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((sigx11_2 * sigx19/sigx3 + sigx8_2 * ewu2/sigx18 - (0.5 *
      sigx9_2 - 0.5 * sigx10_2) * wzdeno * depsi2 * pmusig2/sigx18^2) * prC *
      ewz/sigx3), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((0.5 * (depsi1 * musi1^2/(sigma_sq1)) - depsi1 * muvu1 *
      sigx13_1/sigmastar1) * ewv * sigx13_1/(sigma_sq1) + depsi1 * (mu1/ssq1 -
      (((3 * (mu1) - 2 * (sigx12_1 * (sigma_sq1) * muvu1 * sigmastar1/ssq1^2)) *
        ewv - S * ewu1 * (epsilon)) * sigx12_1 + (0.5 * (ewv/(sigma_sq1)) -
        0.5 * (0.5 * prV1 + ewv/(sigma_sq1))) * prV1 * ewu1 * muvu1/sigmastar1)/ssq1^2)) *
      dmusig1 + (0.5 * (((0.5 * (musi1^2/(sigma_sq1)) - 2) * pmusig1/(sigma_sq1) +
      dmusig1 * sigx13_1) * ewv) + 0.5 * pmusig1) * depsi1 * musi1^2/(sigma_sq1)^2)/(wzdeno *
      pmu1) - ((0.5 * (sigx16_1/(sigma_sq1)) + 0.5 * (((dmusig1 * sigx13_1 -
      wzdeno^2 * pmusig1 * pmu1^2/pmusq1^2) * depsi1 + 0.5 * sigx6_1) * wzdeno *
      pmu1/pmusq1^2)) * ewv + 0.5 * sigx15_1)) * ewz/sqrt(sigma_sq1) + ((((0.5 *
      (depsi2 * musi2^2/(sigma_sq2)) - depsi2 * muvu2 * sigx13_2/sigmastar2) *
      ewv * sigx13_2/(sigma_sq2) + depsi2 * (mu2/ssq2 - (((3 * (mu2) - 2 *
      (sigx12_2 * (sigma_sq2) * muvu2 * sigmastar2/ssq2^2)) * ewv - S * ewu2 *
      (epsilon)) * sigx12_2 + (0.5 * (ewv/(sigma_sq2)) - 0.5 * (0.5 * prV2 +
      ewv/(sigma_sq2))) * prV2 * ewu2 * muvu2/sigmastar2)/ssq2^2)) * dmusig2 +
      (0.5 * (((0.5 * (musi2^2/(sigma_sq2)) - 2) * pmusig2/(sigma_sq2) + dmusig2 *
        sigx13_2) * ewv) + 0.5 * pmusig2) * depsi2 * musi2^2/(sigma_sq2)^2)/pmu2 -
      ((0.5 * (sigx16_2/(sigma_sq2)) + 0.5 * (((dmusig2 * sigx13_2 - pmusig2 *
        pmu2^2/pmusq2^2) * depsi2 + 0.5 * sigx6_2) * pmu2/pmusq2^2)) * ewv +
        0.5 * sigx15_2)) * prC/sqrt(sigma_sq2) - (sigx16_1 * ewz/sqrt(sigma_sq1) +
      sigx16_2 * prC/sqrt(sigma_sq2))^2 * ewv/sigx3) * ewv/sigx3, FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (sigx17 * depsi1 * musi1^2/(sigma_sq1)^2) -
      ((0.5/sqrt(sigma_sq1) - sigx26_1) * ewz + 0.5 * (wzdeno/sqrt(sigma_sq1))) *
        depsi1 * pmu1/pmusq1^2) * pmusig1 + sigx17 * dmusig1 * depsi1 * sigx13_1 -
      ((sigx16_1 * ewz/sqrt(sigma_sq1) + sigx16_2 * prC/sqrt(sigma_sq2)) *
        sigx19/sigx3 + (sigx14_2/(wzdeno * pmu2) - 0.5 * (wzdeno * depsi2 *
        pmusig2 * pmu2/sigx18^2)) * prC/sqrt(sigma_sq2))) * ewv * ewz/sigx3,
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * pmu2 * sqrt(sigma_sq2)) +
      pmu2 * sqrt(sigma_sq2)/sigx18^2) * depsi2 * pmusig2 - (sigx19^2/sigx3 +
      (2 - 2 * (wzdeno * (sigma_sq1) * ewz * pmu1^2/pmusq1^2)) * depsi1 * pmusig1 *
        pmu1 * sqrt(sigma_sq1)/pmusq1^2)) * ewz + sigx17 * depsi1 * pmusig1 -
      prC * depsi2 * pmusig2/sigx18) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# cauchit specification class membership
chessmisfhalfnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/psq2 + ewz1 * depsi1 * pmusig1/psq1)
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx18 <- (pi * sigx2 * ((Wz)^2 + 1))
  sigx19 <- (sigx16_1 * ewz1/sqrt(sigma_sq1) + sigx16_2 * ewz2/sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv) * ewu1/starsq1) *
      ewz1/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 *
      mupsi2/starsq2) * depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv) *
      ewu2/starsq2) * ewz2/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 * depsi1 *
      mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      ewz1/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 *
      depsi2 * mupsi2/sigma_sq2 + (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv +
      1/sigmastar2)/sigma_sq2) * depsi2) * dmusig2) * ewu2/psq2 - (sigx3 *
      sigx11_2/sigx2 + (0.5 * wup2 - 0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) *
      ewz2/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv) *
      sigx13_1/sigma_sq1 - sigx12_1 * depsi1 * ewu1/starsq1^2) * dmusig1 +
      0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
        depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 * (sigx1_1 * pmu1/(sigma_sq1 *
      psq1^2))) * ewz1/sqrt(sigma_sq1) + ((((depsi2 * mupsi2 - depsi2 * musi2/ewv) *
      sigx13_2/sigma_sq2 - sigx12_2 * depsi2 * ewu2/starsq2^2) * dmusig2 +
      0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
        depsi2 * mupsi2/sigma_sq2^2))/pmu2 - 0.5 * (sigx1_2 * pmu2/(sigma_sq2 *
      psq2^2))) * ewz2/sqrt(sigma_sq2) - sigx19 * sigx3/sigx2) * ewv/sigx2,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((sigx1_1/ssqq1 -
    sigx1_2/ssqq2)/sigx18 - pi * sigx3 * ((Wz)^2 + 1) * sigx9/sigx18^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * ewz1/sigx2)) *
      ewz1/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * ewz2 * ewz1/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 * (depsi1 * mupsi1^2/sigma_sq1)) *
      sigx13_1/sigma_sq1 - ((0.5 * (usq1 * ewv/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
      1) * ewv/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu1 *
      sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
      S * (epsilon))) * depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) -
      2) * pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - (0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2) + 0.5 *
      (sigx16_1 * ewu1/sigma_sq1)))/sqrt(sigma_sq1) - sigx19 * sigx11_1/sigx2) *
    ewz1 * ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx11_1 * (1/sigx18 - pi * ((Wz)^2 + 1) * ewz1 * sigx9/sigx18^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((0.5 *
    (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) *
      ewu2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 *
    ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) - 0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) *
    usq2 * ewv * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 *
    sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 +
    sigx8_2 * (0.5 * (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 -
    sigx10_2 * (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 *
    (ewu2/sigma_sq2)) * pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) -
    0.5 * (mu2 * ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
    ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 + 2 *
    ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
    (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * ewz2/sigx2)) *
    ewz2/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 * mupsi2^2/sigma_sq2)) *
      sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
      1) * ewv/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 + mu2 *
      sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
      S * (epsilon))) * depsi2/starsq2^2) * dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) -
      2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) *
      ewu2 + 0.5 * (mu2 * sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - (0.5 * ((sigx10_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) *
      pmu2^2 * sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2) + 0.5 *
      (sigx16_2 * ewu2/sigma_sq2)))/sqrt(sigma_sq2) - sigx19 * sigx11_2/sigx2) *
    ewz2 * ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_2 * (1/sigx18 + pi * ((Wz)^2 + 1) * ewz2 * sigx9/sigx18^2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 *
      sigx13_1/sigmastar1) * ewv * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 -
      (((3 * (mu1) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
        ewv - S * ewu1 * (epsilon)) * sigx12_1 + (0.5 * (ewv/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) *
      dmusig1 + (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 +
      dmusig1 * sigx13_1) * ewv) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      ((0.5 * (sigx16_1/sigma_sq1) + 0.5 * (((dmusig1 * sigx13_1 - pmusig1 *
        pmu1^2/psq1^2) * depsi1 + 0.5 * sigx6_1) * pmu1/psq1^2)) * ewv +
        0.5 * sigx14_1)) * ewz1/sqrt(sigma_sq1) + ((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) -
      depsi2 * musi2 * sigx13_2/sigmastar2) * ewv * sigx13_2/sigma_sq2 + depsi2 *
      (mu2/starsq2 - (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      ((0.5 * (sigx16_2/sigma_sq2) + 0.5 * (((dmusig2 * sigx13_2 - pmusig2 *
        pmu2^2/psq2^2) * depsi2 + 0.5 * sigx6_2) * pmu2/psq2^2)) * ewv +
        0.5 * sigx14_2)) * ewz2/sqrt(sigma_sq2) - sigx19^2 * ewv/sigx2) *
      ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx16_1/sqrt(sigma_sq1) - sigx16_2/sqrt(sigma_sq2))/sigx18 -
      pi * sigx19 * ((Wz)^2 + 1) * sigx9/sigx18^2) * ewv, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx2) + depsi1 * pmusig1/psq1 -
      depsi2 * pmusig2/psq2) * sigx9/sigx18^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# probit specification class membership
chessmisfhalfnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/psq2 + depsi1 * pmusig1 * pwZ/psq1)
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx19 <- (sigx16_1 * pwZ/sqrt(sigma_sq1) + sigx16_2 * (1 - pwZ)/sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv) * ewu1/starsq1) *
      pwZ/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
      depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv) * ewu2/starsq2) *
      (1 - pwZ)/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 * depsi1 *
      mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      pwZ/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 *
      depsi2 * mupsi2/sigma_sq2 + (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv +
      1/sigmastar2)/sigma_sq2) * depsi2) * dmusig2) * ewu2/psq2 - (sigx3 *
      sigx11_2/sigx2 + (0.5 * wup2 - 0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) *
      (1 - pwZ)/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv) *
      sigx13_1/sigma_sq1 - sigx12_1 * depsi1 * ewu1/starsq1^2) * dmusig1 +
      0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
        depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 * (sigx1_1 * pmu1/(sigma_sq1 *
      psq1^2))) * pwZ/sqrt(sigma_sq1) + ((((depsi2 * mupsi2 - depsi2 * musi2/ewv) *
      sigx13_2/sigma_sq2 - sigx12_2 * depsi2 * ewu2/starsq2^2) * dmusig2 +
      0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
        depsi2 * mupsi2/sigma_sq2^2))/pmu2 - 0.5 * (sigx1_2 * pmu2/(sigma_sq2 *
      psq2^2))) * (1 - pwZ)/sqrt(sigma_sq2) - sigx19 * sigx3/sigx2) * ewv/sigx2,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx1_1/ssqq1 -
    (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * pwZ/sigx2)) *
      pwZ/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * (1 - pwZ) * pwZ/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 * (depsi1 * mupsi1^2/sigma_sq1)) *
      sigx13_1/sigma_sq1 - ((0.5 * (usq1 * ewv/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
      1) * ewv/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu1 *
      sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
      S * (epsilon))) * depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) -
      2) * pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - (0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2) + 0.5 *
      (sigx16_1 * ewu1/sigma_sq1)))/sqrt(sigma_sq1) - sigx19 * sigx11_1/sigx2) *
    pwZ * ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx11_1 * (1 - sigx9 * pwZ/sigx2) * dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((0.5 *
    (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) *
      ewu2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 *
    ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) - 0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) *
    usq2 * ewv * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 *
    sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 +
    sigx8_2 * (0.5 * (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 -
    sigx10_2 * (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 *
    (ewu2/sigma_sq2)) * pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) -
    0.5 * (mu2 * ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
    ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 + 2 *
    ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
    (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * (1 - pwZ)/sigx2)) *
    (1 - pwZ)/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 * mupsi2^2/sigma_sq2)) *
      sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
      1) * ewv/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 + mu2 *
      sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
      S * (epsilon))) * depsi2/starsq2^2) * dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) -
      2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) *
      ewu2 + 0.5 * (mu2 * sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - (0.5 * ((sigx10_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) *
      pmu2^2 * sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2) + 0.5 *
      (sigx16_2 * ewu2/sigma_sq2)))/sqrt(sigma_sq2) - sigx19 * sigx11_2/sigx2) *
    (1 - pwZ) * ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_2 * ((1 - pwZ) * sigx9/sigx2 + 1) * dwZ/sigx2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 *
      sigx13_1/sigmastar1) * ewv * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 -
      (((3 * (mu1) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
        ewv - S * ewu1 * (epsilon)) * sigx12_1 + (0.5 * (ewv/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) *
      dmusig1 + (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 +
      dmusig1 * sigx13_1) * ewv) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      ((0.5 * (sigx16_1/sigma_sq1) + 0.5 * (((dmusig1 * sigx13_1 - pmusig1 *
        pmu1^2/psq1^2) * depsi1 + 0.5 * sigx6_1) * pmu1/psq1^2)) * ewv +
        0.5 * sigx14_1)) * pwZ/sqrt(sigma_sq1) + ((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) -
      depsi2 * musi2 * sigx13_2/sigmastar2) * ewv * sigx13_2/sigma_sq2 + depsi2 *
      (mu2/starsq2 - (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      ((0.5 * (sigx16_2/sigma_sq2) + 0.5 * (((dmusig2 * sigx13_2 - pmusig2 *
        pmu2^2/psq2^2) * depsi2 + 0.5 * sigx6_2) * pmu2/psq2^2)) * ewv +
        0.5 * sigx14_2)) * (1 - pwZ)/sqrt(sigma_sq2) - sigx19^2 * ewv/sigx2) *
      ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx16_1/sqrt(sigma_sq1) - (sigx19 * sigx9/sigx2 +
      sigx16_2/sqrt(sigma_sq2))) * dwZ * ewv/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx9 * dwZ/sigx2 + Wz) * sigx9 * dwZ/sigx2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# cloglog specification class membership
chessmisfhalfnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  mu1 <- 0
  mu2 <- 0
  musi1 <- (mu1 * ewv - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/psq1 + depsi2 * prZ * pmusig2/psq2)
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx19 <- (sigx16_1 * (1 - prZ)/sqrt(sigma_sq1) + sigx16_2 * prZ/sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv) * ewu1/starsq1) *
      (1 - prZ)/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 *
      mupsi2/starsq2) * depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv) *
      ewu2/starsq2) * prZ/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 *
      ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 * depsi1 *
      mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      (1 - prZ)/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 +
      dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 *
      depsi2 * mupsi2/sigma_sq2 + (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv +
      1/sigmastar2)/sigma_sq2) * depsi2) * dmusig2) * ewu2/psq2 - (sigx3 *
      sigx11_2/sigx2 + (0.5 * wup2 - 0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) *
      prZ/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv) *
      sigx13_1/sigma_sq1 - sigx12_1 * depsi1 * ewu1/starsq1^2) * dmusig1 +
      0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
        depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 * (sigx1_1 * pmu1/(sigma_sq1 *
      psq1^2))) * (1 - prZ)/sqrt(sigma_sq1) + ((((depsi2 * mupsi2 - depsi2 *
      musi2/ewv) * sigx13_2/sigma_sq2 - sigx12_2 * depsi2 * ewu2/starsq2^2) *
      dmusig2 + 0.5 * (((mupsi2^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu2 *
      mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2))/pmu2 - 0.5 * (sigx1_2 *
      pmu2/(sigma_sq2 * psq2^2))) * prZ/sqrt(sigma_sq2) - sigx19 * sigx3/sigx2) *
      ewv/sigx2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (sigx1_1/ssqq1 -
    (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * (1 - prZ)/sigx2)) *
      (1 - prZ)/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * prZ * (1 - prZ)/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 * (depsi1 * mupsi1^2/sigma_sq1)) *
      sigx13_1/sigma_sq1 - ((0.5 * (usq1 * ewv/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
      1) * ewv/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu1 *
      sigx7_1 - sigx12_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
      S * (epsilon))) * depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) -
      2) * pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - (0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2) + 0.5 *
      (sigx16_1 * ewu1/sigma_sq1)))/sqrt(sigma_sq1) - sigx19 * sigx11_1/sigx2) *
    (1 - prZ) * ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx11_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((((0.5 *
    (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) *
      ewu2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 *
    ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) - 0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) *
    usq2 * ewv * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 *
    sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 +
    sigx8_2 * (0.5 * (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 -
    sigx10_2 * (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 *
    (ewu2/sigma_sq2)) * pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) -
    0.5 * (mu2 * ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
    ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 + 2 *
    ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
    (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * prZ/sigx2)) * prZ/sigx2,
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 * mupsi2^2/sigma_sq2)) *
      sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
      1) * ewv/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 + mu2 *
      sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
      S * (epsilon))) * depsi2/starsq2^2) * dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) -
      2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) *
      ewu2 + 0.5 * (mu2 * sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - (0.5 * ((sigx10_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) *
      pmu2^2 * sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2) + 0.5 *
      (sigx16_2 * ewu2/sigma_sq2)))/sqrt(sigma_sq2) - sigx19 * sigx11_2/sigx2) *
    prZ * ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_2 * (sigx9 * prZ/sigx2 + 1) * prZ * ewz/sigx2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 *
      sigx13_1/sigmastar1) * ewv * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 -
      (((3 * (mu1) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
        ewv - S * ewu1 * (epsilon)) * sigx12_1 + (0.5 * (ewv/sigma_sq1) -
        0.5 * (0.5 * vsq1 + ewv/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) *
      dmusig1 + (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 +
      dmusig1 * sigx13_1) * ewv) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      ((0.5 * (sigx16_1/sigma_sq1) + 0.5 * (((dmusig1 * sigx13_1 - pmusig1 *
        pmu1^2/psq1^2) * depsi1 + 0.5 * sigx6_1) * pmu1/psq1^2)) * ewv +
        0.5 * sigx14_1)) * (1 - prZ)/sqrt(sigma_sq1) + ((((0.5 * (depsi2 *
      mupsi2^2/sigma_sq2) - depsi2 * musi2 * sigx13_2/sigmastar2) * ewv * sigx13_2/sigma_sq2 +
      depsi2 * (mu2/starsq2 - (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 *
        sigmastar2/starsq2^2)) * ewv - S * ewu2 * (epsilon)) * sigx12_2 +
        (0.5 * (ewv/sigma_sq2) - 0.5 * (0.5 * vsq2 + ewv/sigma_sq2)) * vsq2 *
          ewu2 * musi2/sigmastar2)/starsq2^2)) * dmusig2 + (0.5 * (((0.5 *
      (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 + dmusig2 * sigx13_2) *
      ewv) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 - ((0.5 *
      (sigx16_2/sigma_sq2) + 0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/psq2^2) *
      depsi2 + 0.5 * sigx6_2) * pmu2/psq2^2)) * ewv + 0.5 * sigx14_2)) * prZ/sqrt(sigma_sq2) -
      sigx19^2 * ewv/sigx2) * ewv/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx16_1/sqrt(sigma_sq1) - (sigx19 * sigx9/sigx2 +
      sigx16_2/sqrt(sigma_sq2))) * prZ * ewv * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx9 * prZ/sigx2 + 1) * ewz) * sigx9 *
      prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for misf halfnormal-normal distribution
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
misfhalfnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfhalfnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfhalfnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfhalfnormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfhalfnormlike_logit,
    grad = cgradmisfhalfnormlike_logit, hess = chessmisfhalfnormlike_logit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfhalfnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfhalfnormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfhalfnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfhalfnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfhalfnormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfhalfnormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfhalfnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfhalfnormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfhalfnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

# cauchit specification class membership
misfhalfnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfhalfnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfhalfnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfhalfnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfhalfnormlike_cauchit,
    grad = cgradmisfhalfnormlike_cauchit, hess = chessmisfhalfnormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfhalfnormlike_cauchit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfhalfnormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfhalfnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfhalfnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfhalfnormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfhalfnormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfhalfnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfhalfnormlike_cauchit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfhalfnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

# probit specification class membership
misfhalfnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfhalfnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfhalfnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfhalfnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfhalfnormlike_probit,
    grad = cgradmisfhalfnormlike_probit, hess = chessmisfhalfnormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfhalfnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfhalfnormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfhalfnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfhalfnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfhalfnormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfhalfnormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfhalfnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfhalfnormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfhalfnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

# cloglog specification class membership
misfhalfnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfhalfnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfhalfnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfhalfnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfhalfnormlike_cloglog,
    grad = cgradmisfhalfnormlike_cloglog, hess = chessmisfhalfnormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfhalfnormlike_cloglog(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfhalfnormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfhalfnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfhalfnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfhalfnormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfhalfnormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfhalfnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfhalfnormlike_cloglog(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfhalfnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initHalf = initHalf))
}

# Conditional efficiencies estimation ----------
#' efficiencies for misf halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# logit specification class membership
cmisfhalfnormeff_logit <- function(object, level) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

# cauchit specification class membership
cmisfhalfnormeff_cauchit <- function(object, level) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

# probit specification class membership
cmisfhalfnormeff_probit <- function(object, level) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

# cloglog specification class membership
cmisfhalfnormeff_cloglog <- function(object, level) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }

  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for misf halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
# logit specification class membership
cmisfmarghalfnorm_Eu_logit <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarghalfnorm_Vu_logit <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cauchit specification class membership
cmisfmarghalfnorm_Eu_cauchit <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarghalfnorm_Vu_cauchit <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# probit specification class membership
cmisfmarghalfnorm_Eu_probit <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarghalfnorm_Vu_probit <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cloglog specification class membership
cmisfmarghalfnorm_Eu_cloglog <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarghalfnorm_Vu_cloglog <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv))) * pnorm(mustar2/sigmastar2)
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
