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
# Convolution: rayleigh - normal                                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for misf rayleigh-normal distribution
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
cmisfraynormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# cauchit specification class membership
cmisfraynormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# probit specification class membership
cmisfraynormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# cloglog specification class membership
cmisfraynormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for misf rayleigh-normal distribution
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
cstmisfraynorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart,
  initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- cstraynorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initRay <- NULL

  } else {
    cat("Initialization: SFA + rayleigh - normal distributions...\n")
    initRay <- maxLik::maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradraynormlike, method = initAlg,
      control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], vHvar = vHvar[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar)
    Esti <- initRay$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("MI_", colnames(Zvar)))
  return(list(StartVal = StartVal, initRay = initRay))
}

# Gradient of the likelihood function ----------
#' gradient for misf rayleigh-normal distribution
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
cgradmisfraynormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
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
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  suvq1 <- (sigmastar1/ewv - ewu1/ssq1)
  suvq2 <- (sigmastar2/ewv - ewu2/ssq2)
  sigx2_1 <- (pmusig1 + S * dmusig1 * suvq1 * (epsilon))
  sigx2_2 <- (pmusig2 + S * dmusig2 * suvq2 * (epsilon))
  wusq1 <- (1 - ewu1/(sigma_sq1))
  wusq2 <- (1 - ewu2/(sigma_sq2))
  wvsq1 <- (1 - ewv/(sigma_sq1))
  wvsq2 <- (1 - ewv/(sigma_sq2))
  sigx3_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/(sigma_sq1))
  sigx3_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/(sigma_sq2))
  sigx4_1 <- (ewu1 * sigx2_1/(sigma_sq1) + S * wusq1 * sigx3_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx2_2/(sigma_sq2) + S * wusq2 * sigx3_2 *
    (epsilon)/ewv)
  sigx5 <- (prC * sigx1_2 * sigx4_2 * sigmastar2/ewu2 + sigx1_1 *
    sigx4_1 * ewz * sigmastar1/(wzdeno * ewu1))
  sigx6 <- (prC * sigx3_2 * sigx1_2 * sigmastar2/ewu2 + sigx3_1 *
    sigx1_1 * ewz * sigmastar1/(wzdeno * ewu1))
  dvsq1 <- (dmusig1 * ewv/sigmastar1)
  dvsq2 <- (dmusig2 * ewv/sigmastar2)
  sigx7_1 <- (0.5 * dvsq1 - S * pmusig1 * (epsilon))
  sigx7_2 <- (0.5 * dvsq2 - S * pmusig2 * (epsilon))
  wuvsq1 <- (wusq1 * ewv/sigmastar1)
  wuvsq2 <- (wusq2 * ewv/sigmastar2)
  s8sig1 <- (0.5 * wuvsq1 + sigmastar1)
  s8sig2 <- (0.5 * wuvsq2 + sigmastar2)
  sigx8_1 <- (1/ssq1 - s8sig1 * ewu1/ssq1^2)
  sigx8_2 <- (1/ssq2 - s8sig2 * ewu2/ssq2^2)
  sigx9_1 <- (sigx7_1 * wusq1 + S^2 * sigx8_1 * sigx3_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx9_2 <- (sigx7_2 * wusq2 + S^2 * sigx8_2 * sigx3_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx10_1 <- (wusq1 * sigx3_1 * ewv/sigmastar1)
  sigx10_2 <- (wusq2 * sigx3_2 * ewv/sigmastar2)
  sigx11_1 <- (sigx9_1 * sigmastar1 + 0.5 * sigx10_1)
  sigx11_2 <- (sigx9_2 * sigmastar2 + 0.5 * sigx10_2)
  sigx12_1 <- (wzdeno * (sigma_sq1))
  sigx13_1 <- (sigx11_1/sigx12_1 - wzdeno * sigx3_1 * ewu1 *
    sigmastar1/(wzdeno * ewu1)^2)
  sigx13_2 <- (sigx11_2/(sigma_sq2) - sigx3_2 * sigmastar2/ewu2)
  sigx14_1 <- (wvsq1 * dmusig1/sigmastar1)
  sigx14_2 <- (wvsq2 * dmusig2/sigmastar2)
  sigx15_1 <- (0.5 * sigx14_1 + S * pmusig1 * (epsilon)/(sigma_sq1))
  sigx15_2 <- (0.5 * sigx14_2 + S * pmusig2 * (epsilon)/(sigma_sq2))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx17_1 <- (wvsq1 * ewu1/sigmastar1)
  sigx17_2 <- (wvsq2 * ewu2/sigmastar2)
  sigx18_1 <- (0.5 * sigx17_1 + sigmastar1)
  sigx18_2 <- (0.5 * sigx17_2 + sigmastar2)
  sigx19_1 <- (ssq1^2 * (sigma_sq1) * sigmastar1)
  sigx19_2 <- (ssq2^2 * (sigma_sq2) * sigmastar2)
  sigx20_1 <- (2 * sigx16 - S^2 * sigx18_1 * ewu1^2 * (epsilon)^2/sigx19_1)
  sigx20_2 <- (2 * sigx16 - S^2 * sigx18_2 * ewu2^2 * (epsilon)^2/sigx19_2)
  sigx21_1 <- (sigx15_1 * ewu1/(sigma_sq1) + sigx20_1 * sigx3_1)
  sigx21_2 <- (sigx15_2 * ewu2/(sigma_sq2) + sigx20_2 * sigx3_2)
  sigx22_1 <- (wvsq1 * sigx3_1 * ewu1/ssq1)
  sigx22_2 <- (wvsq2 * sigx3_2 * ewu2/ssq2)
  sigx23 <- (1/(wzdeno * ewu1) - ewu1 * ewz/(wzdeno * ewu1)^2)
  sigx24_1 <- (sigx21_1 * sigmastar1 + 0.5 * sigx22_1)
  sigx24_2 <- (sigx21_2 * sigmastar2 + 0.5 * sigx22_2)
  sigx25 <- (sigx24_1 * sigx1_1 * ewz/(wzdeno * ewu1) + sigx24_2 *
    prC * sigx1_2/ewu2)
  sigx26 <- (sigx23 * sigx3_1 * sigx1_1 * sigmastar1 - prC *
    sigx3_2 * sigx1_2 * sigmastar2/(wzdeno * ewu2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx6,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13_1 *
    sigx1_1 * ewz/sigx6, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = sigx13_2 * prC * sigx1_2/sigx6, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = (sigx25 * ewv/sigx6 - 0.5), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx26 * ewz/sigx6, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cauchit specification class membership
cgradmisfraynormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  musig1 <- (S * ewu1 * (epsilon)/starsq1)
  musig2 <- (S * ewu2 * (epsilon)/starsq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  wuv1 <- (sigmastar1/ewv - ewu1/starsq1)
  wuv2 <- (sigmastar2/ewv - ewu2/starsq2)
  dvs1 <- (dmusig1 * ewv/sigmastar1)
  dvs2 <- (dmusig2 * ewv/sigmastar2)
  dpsv1 <- (0.5 * dvs1 - S * pmusig1 * (epsilon))
  dpsv2 <- (0.5 * dvs2 - S * pmusig2 * (epsilon))
  vepsi <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx2_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/sigma_sq1)
  sigx2_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/sigma_sq2)
  sigx3_1 <- (pmusig1 + S * dmusig1 * wuv1 * (epsilon))
  sigx3_2 <- (pmusig2 + S * dmusig2 * wuv2 * (epsilon))
  sigx4_1 <- (ewu1 * sigx3_1/sigma_sq1 + S * usq1 * sigx2_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx3_2/sigma_sq2 + S * usq2 * sigx2_2 *
    (epsilon)/ewv)
  sigx5_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx5_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx6_1 <- (1/starsq1 - sigx5_1 * ewu1/starsq1^2)
  sigx6_2 <- (1/starsq2 - sigx5_2 * ewu2/starsq2^2)
  sigx7_1 <- (dpsv1 * usq1 + S^2 * sigx6_1 * sigx2_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx7_2 <- (dpsv2 * usq2 + S^2 * sigx6_2 * sigx2_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx8_1 <- (usq1 * sigx2_1 * ewv/sigmastar1)
  sigx8_2 <- (usq2 * sigx2_2 * ewv/sigmastar2)
  sigx9_1 <- (sigx7_1 * sigmastar1 + 0.5 * sigx8_1)
  sigx9_2 <- (sigx7_2 * sigmastar2 + 0.5 * sigx8_2)
  sigx10_1 <- (vsq1 * dmusig1/sigmastar1)
  sigx10_2 <- (vsq2 * dmusig2/sigmastar2)
  sigx11_1 <- (0.5 * sigx10_1 + S * pmusig1 * (epsilon)/sigma_sq1)
  sigx11_2 <- (0.5 * sigx10_2 + S * pmusig2 * (epsilon)/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (starsq1^2 * sigma_sq1 * sigmastar1)
  sigx13_2 <- (starsq2^2 * sigma_sq2 * sigmastar2)
  sigx14_1 <- (2 * vepsi - S^2 * sigx12_1 * ewu1^2 * (epsilon)^2/sigx13_1)
  sigx14_2 <- (2 * vepsi - S^2 * sigx12_2 * ewu2^2 * (epsilon)^2/sigx13_2)
  sigx15_1 <- (sigx11_1 * ewu1/sigma_sq1 + sigx14_1 * sigx2_1)
  sigx15_2 <- (sigx11_2 * ewu2/sigma_sq2 + sigx14_2 * sigx2_2)
  sigx16_1 <- (sigx15_1 * sigmastar1 + 0.5 * (vsq1 * sigx2_1 *
    ewu1/starsq1))
  sigx16_2 <- (sigx15_2 * sigmastar2 + 0.5 * (vsq2 * sigx2_2 *
    ewu2/starsq2))
  sigx17_1 <- (sigx9_1/sigma_sq1 - sigx2_1 * sigmastar1/ewu1)
  sigx17_2 <- (sigx9_2/sigma_sq2 - sigx2_2 * sigmastar2/ewu2)
  sigx18 <- (ewz2 * sigx1_2 * sigx4_2 * sigmastar2/ewu2 + ewz1 *
    sigx1_1 * sigx4_1 * sigmastar1/ewu1)
  sigx19 <- (ewz2 * sigx2_2 * sigx1_2 * sigmastar2/ewu2 + ewz1 *
    sigx2_1 * sigx1_1 * sigmastar1/ewu1)
  sigx20 <- (sigx16_1 * ewz1 * sigx1_1/ewu1 + sigx16_2 * ewz2 *
    sigx1_2/ewu2)
  sigx21 <- (sigx2_1 * sigx1_1 * sigmastar1/ewu1 - sigx2_2 *
    sigx1_2 * sigmastar2/ewu2)
  sigx22 <- (pi * sigx19 * ((Wz)^2 + 1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx18/sigx19,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx17_1 *
    ewz1 * sigx1_1/sigx19, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = sigx17_2 * ewz2 * sigx1_2/sigx19, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (sigx20 * ewv/sigx19 -
      0.5), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx21/sigx22,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# probit specification class membership
cgradmisfraynormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  musig1 <- (S * ewu1 * (epsilon)/starsq1)
  musig2 <- (S * ewu2 * (epsilon)/starsq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  wuv1 <- (sigmastar1/ewv - ewu1/starsq1)
  wuv2 <- (sigmastar2/ewv - ewu2/starsq2)
  dvs1 <- (dmusig1 * ewv/sigmastar1)
  dvs2 <- (dmusig2 * ewv/sigmastar2)
  dpsv1 <- (0.5 * dvs1 - S * pmusig1 * (epsilon))
  dpsv2 <- (0.5 * dvs2 - S * pmusig2 * (epsilon))
  vepsi <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx2_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/sigma_sq1)
  sigx2_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/sigma_sq2)
  sigx3_1 <- (pmusig1 + S * dmusig1 * wuv1 * (epsilon))
  sigx3_2 <- (pmusig2 + S * dmusig2 * wuv2 * (epsilon))
  sigx4_1 <- (ewu1 * sigx3_1/sigma_sq1 + S * usq1 * sigx2_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx3_2/sigma_sq2 + S * usq2 * sigx2_2 *
    (epsilon)/ewv)
  sigx5_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx5_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx6_1 <- (1/starsq1 - sigx5_1 * ewu1/starsq1^2)
  sigx6_2 <- (1/starsq2 - sigx5_2 * ewu2/starsq2^2)
  sigx7_1 <- (dpsv1 * usq1 + S^2 * sigx6_1 * sigx2_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx7_2 <- (dpsv2 * usq2 + S^2 * sigx6_2 * sigx2_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx8_1 <- (usq1 * sigx2_1 * ewv/sigmastar1)
  sigx8_2 <- (usq2 * sigx2_2 * ewv/sigmastar2)
  sigx9_1 <- (sigx7_1 * sigmastar1 + 0.5 * sigx8_1)
  sigx9_2 <- (sigx7_2 * sigmastar2 + 0.5 * sigx8_2)
  sigx10_1 <- (vsq1 * dmusig1/sigmastar1)
  sigx10_2 <- (vsq2 * dmusig2/sigmastar2)
  sigx11_1 <- (0.5 * sigx10_1 + S * pmusig1 * (epsilon)/sigma_sq1)
  sigx11_2 <- (0.5 * sigx10_2 + S * pmusig2 * (epsilon)/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (starsq1^2 * sigma_sq1 * sigmastar1)
  sigx13_2 <- (starsq2^2 * sigma_sq2 * sigmastar2)
  sigx14_1 <- (2 * vepsi - S^2 * sigx12_1 * ewu1^2 * (epsilon)^2/sigx13_1)
  sigx14_2 <- (2 * vepsi - S^2 * sigx12_2 * ewu2^2 * (epsilon)^2/sigx13_2)
  sigx15_1 <- (sigx11_1 * ewu1/sigma_sq1 + sigx14_1 * sigx2_1)
  sigx15_2 <- (sigx11_2 * ewu2/sigma_sq2 + sigx14_2 * sigx2_2)
  sigx16_1 <- (sigx15_1 * sigmastar1 + 0.5 * (vsq1 * sigx2_1 *
    ewu1/starsq1))
  sigx16_2 <- (sigx15_2 * sigmastar2 + 0.5 * (vsq2 * sigx2_2 *
    ewu2/starsq2))
  sigx17_1 <- (sigx9_1/sigma_sq1 - sigx2_1 * sigmastar1/ewu1)
  sigx17_2 <- (sigx9_2/sigma_sq2 - sigx2_2 * sigmastar2/ewu2)
  sigx18 <- ((1 - pwZ) * sigx1_2 * sigx4_2 * sigmastar2/ewu2 +
    sigx1_1 * sigx4_1 * pwZ * sigmastar1/ewu1)
  sigx19 <- ((1 - pwZ) * sigx2_2 * sigx1_2 * sigmastar2/ewu2 +
    sigx2_1 * sigx1_1 * pwZ * sigmastar1/ewu1)
  sigx20 <- (sigx16_1 * sigx1_1 * pwZ/ewu1 + sigx16_2 * (1 -
    pwZ) * sigx1_2/ewu2)
  sigx21 <- (sigx2_1 * sigx1_1 * sigmastar1/ewu1 - sigx2_2 *
    sigx1_2 * sigmastar2/ewu2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx18/sigx19,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx17_1 *
    pwZ * sigx1_1/sigx19, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = sigx17_2 * (1 - pwZ) * sigx1_2/sigx19, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (sigx20 * ewv/sigx19 -
      0.5), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx21 *
      dwZ/sigx19, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}


# cloglog specification class membership
cgradmisfraynormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  musig1 <- (S * ewu1 * (epsilon)/starsq1)
  musig2 <- (S * ewu2 * (epsilon)/starsq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  wuv1 <- (sigmastar1/ewv - ewu1/starsq1)
  wuv2 <- (sigmastar2/ewv - ewu2/starsq2)
  dvs1 <- (dmusig1 * ewv/sigmastar1)
  dvs2 <- (dmusig2 * ewv/sigmastar2)
  dpsv1 <- (0.5 * dvs1 - S * pmusig1 * (epsilon))
  dpsv2 <- (0.5 * dvs2 - S * pmusig2 * (epsilon))
  vepsi <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx2_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/sigma_sq1)
  sigx2_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/sigma_sq2)
  sigx3_1 <- (pmusig1 + S * dmusig1 * wuv1 * (epsilon))
  sigx3_2 <- (pmusig2 + S * dmusig2 * wuv2 * (epsilon))
  sigx4_1 <- (ewu1 * sigx3_1/sigma_sq1 + S * usq1 * sigx2_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx3_2/sigma_sq2 + S * usq2 * sigx2_2 *
    (epsilon)/ewv)
  sigx5_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx5_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx6_1 <- (1/starsq1 - sigx5_1 * ewu1/starsq1^2)
  sigx6_2 <- (1/starsq2 - sigx5_2 * ewu2/starsq2^2)
  sigx7_1 <- (dpsv1 * usq1 + S^2 * sigx6_1 * sigx2_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx7_2 <- (dpsv2 * usq2 + S^2 * sigx6_2 * sigx2_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx8_1 <- (usq1 * sigx2_1 * ewv/sigmastar1)
  sigx8_2 <- (usq2 * sigx2_2 * ewv/sigmastar2)
  sigx9_1 <- (sigx7_1 * sigmastar1 + 0.5 * sigx8_1)
  sigx9_2 <- (sigx7_2 * sigmastar2 + 0.5 * sigx8_2)
  sigx10_1 <- (vsq1 * dmusig1/sigmastar1)
  sigx10_2 <- (vsq2 * dmusig2/sigmastar2)
  sigx11_1 <- (0.5 * sigx10_1 + S * pmusig1 * (epsilon)/sigma_sq1)
  sigx11_2 <- (0.5 * sigx10_2 + S * pmusig2 * (epsilon)/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (starsq1^2 * sigma_sq1 * sigmastar1)
  sigx13_2 <- (starsq2^2 * sigma_sq2 * sigmastar2)
  sigx14_1 <- (2 * vepsi - S^2 * sigx12_1 * ewu1^2 * (epsilon)^2/sigx13_1)
  sigx14_2 <- (2 * vepsi - S^2 * sigx12_2 * ewu2^2 * (epsilon)^2/sigx13_2)
  sigx15_1 <- (sigx11_1 * ewu1/sigma_sq1 + sigx14_1 * sigx2_1)
  sigx15_2 <- (sigx11_2 * ewu2/sigma_sq2 + sigx14_2 * sigx2_2)
  sigx16_1 <- (sigx15_1 * sigmastar1 + 0.5 * (vsq1 * sigx2_1 *
    ewu1/starsq1))
  sigx16_2 <- (sigx15_2 * sigmastar2 + 0.5 * (vsq2 * sigx2_2 *
    ewu2/starsq2))
  sigx17_1 <- (sigx9_1/sigma_sq1 - sigx2_1 * sigmastar1/ewu1)
  sigx17_2 <- (sigx9_2/sigma_sq2 - sigx2_2 * sigmastar2/ewu2)
  sigx18 <- ((1 - prZ) * sigx1_1 * sigx4_1 * sigmastar1/ewu1 +
    prZ * sigx1_2 * sigx4_2 * sigmastar2/ewu2)
  sigx19 <- ((1 - prZ) * sigx2_1 * sigx1_1 * sigmastar1/ewu1 +
    sigx2_2 * prZ * sigx1_2 * sigmastar2/ewu2)
  sigx20 <- (sigx16_1 * (1 - prZ) * sigx1_1/ewu1 + sigx16_2 *
    prZ * sigx1_2/ewu2)
  sigx21 <- (sigx2_1 * sigx1_1 * sigmastar1/ewu1 - sigx2_2 *
    sigx1_2 * sigmastar2/ewu2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx18/sigx19,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx17_1 *
    (1 - prZ) * sigx1_1/sigx19, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = sigx17_2 * prZ * sigx1_2/sigx19,
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx20 *
    ewv/sigx19 - 0.5), FUN = "*"), sweep(Zvar, MARGIN = 1,
    STATS = sigx21 * prZ * ewz/sigx19, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for misf rayleigh-normal distribution
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
chessmisfraynormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
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
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  suvq1 <- (sigmastar1/ewv - ewu1/ssq1)
  suvq2 <- (sigmastar2/ewv - ewu2/ssq2)
  sigx2_1 <- (pmusig1 + S * dmusig1 * suvq1 * (epsilon))
  sigx2_2 <- (pmusig2 + S * dmusig2 * suvq2 * (epsilon))
  wusq1 <- (1 - ewu1/(sigma_sq1))
  wusq2 <- (1 - ewu2/(sigma_sq2))
  wvsq1 <- (1 - ewv/(sigma_sq1))
  wvsq2 <- (1 - ewv/(sigma_sq2))
  sigx3_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/(sigma_sq1))
  sigx3_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/(sigma_sq2))
  sigx4_1 <- (ewu1 * sigx2_1/(sigma_sq1) + S * wusq1 * sigx3_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx2_2/(sigma_sq2) + S * wusq2 * sigx3_2 *
    (epsilon)/ewv)
  sigx5 <- (prC * sigx1_2 * sigx4_2 * sigmastar2/ewu2 + sigx1_1 *
    sigx4_1 * ewz * sigmastar1/(wzdeno * ewu1))
  sigx6 <- (prC * sigx3_2 * sigx1_2 * sigmastar2/ewu2 + sigx3_1 *
    sigx1_1 * ewz * sigmastar1/(wzdeno * ewu1))
  dvsq1 <- (dmusig1 * ewv/sigmastar1)
  dvsq2 <- (dmusig2 * ewv/sigmastar2)
  sigx7_1 <- (0.5 * dvsq1 - S * pmusig1 * (epsilon))
  sigx7_2 <- (0.5 * dvsq2 - S * pmusig2 * (epsilon))
  wuvsq1 <- (wusq1 * ewv/sigmastar1)
  wuvsq2 <- (wusq2 * ewv/sigmastar2)
  s8sig1 <- (0.5 * wuvsq1 + sigmastar1)
  s8sig2 <- (0.5 * wuvsq2 + sigmastar2)
  sigx8_1 <- (1/ssq1 - s8sig1 * ewu1/ssq1^2)
  sigx8_2 <- (1/ssq2 - s8sig2 * ewu2/ssq2^2)
  sigx9_1 <- (sigx7_1 * wusq1 + S^2 * sigx8_1 * sigx3_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx9_2 <- (sigx7_2 * wusq2 + S^2 * sigx8_2 * sigx3_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx10_1 <- (wusq1 * sigx3_1 * ewv/sigmastar1)
  sigx10_2 <- (wusq2 * sigx3_2 * ewv/sigmastar2)
  sigx11_1 <- (sigx9_1 * sigmastar1 + 0.5 * sigx10_1)
  sigx11_2 <- (sigx9_2 * sigmastar2 + 0.5 * sigx10_2)
  sigx12_1 <- (wzdeno * (sigma_sq1))
  sigx13_1 <- (sigx11_1/sigx12_1 - wzdeno * sigx3_1 * ewu1 *
    sigmastar1/(wzdeno * ewu1)^2)
  sigx13_2 <- (sigx11_2/(sigma_sq2) - sigx3_2 * sigmastar2/ewu2)
  sigx14_1 <- (wvsq1 * dmusig1/sigmastar1)
  sigx14_2 <- (wvsq2 * dmusig2/sigmastar2)
  sigx15_1 <- (0.5 * sigx14_1 + S * pmusig1 * (epsilon)/(sigma_sq1))
  sigx15_2 <- (0.5 * sigx14_2 + S * pmusig2 * (epsilon)/(sigma_sq2))
  sigx16 <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx17_1 <- (wvsq1 * ewu1/sigmastar1)
  sigx17_2 <- (wvsq2 * ewu2/sigmastar2)
  sigx18_1 <- (0.5 * sigx17_1 + sigmastar1)
  sigx18_2 <- (0.5 * sigx17_2 + sigmastar2)
  sigx19_1 <- (ssq1^2 * (sigma_sq1) * sigmastar1)
  sigx19_2 <- (ssq2^2 * (sigma_sq2) * sigmastar2)
  sigx20_1 <- (2 * sigx16 - S^2 * sigx18_1 * ewu1^2 * (epsilon)^2/sigx19_1)
  sigx20_2 <- (2 * sigx16 - S^2 * sigx18_2 * ewu2^2 * (epsilon)^2/sigx19_2)
  sigx21_1 <- (sigx15_1 * ewu1/(sigma_sq1) + sigx20_1 * sigx3_1)
  sigx21_2 <- (sigx15_2 * ewu2/(sigma_sq2) + sigx20_2 * sigx3_2)
  sigx22_1 <- (wvsq1 * sigx3_1 * ewu1/ssq1)
  sigx22_2 <- (wvsq2 * sigx3_2 * ewu2/ssq2)
  sigx23 <- (1/(wzdeno * ewu1) - ewu1 * ewz/(wzdeno * ewu1)^2)
  sigx24_1 <- (sigx21_1 * sigmastar1 + 0.5 * sigx22_1)
  sigx24_2 <- (sigx21_2 * sigmastar2 + 0.5 * sigx22_2)
  sigx25 <- (sigx24_1 * sigx1_1 * ewz/(wzdeno * ewu1) + sigx24_2 *
    prC * sigx1_2/ewu2)
  sigx26 <- (sigx23 * sigx3_1 * sigx1_1 * sigmastar1 - prC *
    sigx3_2 * sigx1_2 * sigmastar2/(wzdeno * ewu2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((S^2 * ewu1 * (epsilon)^2/((sigma_sq1) *
      ewv) - 1) * suvq1 + ewu1/ssq1) * dmusig1 * ewu1/(sigma_sq1) +
      wusq1 * (S * ((2 * (S * dmusig1 * suvq1 * (epsilon)) +
        3 * pmusig1) * ewu1/(sigma_sq1) + S * wusq1 *
        sigx3_1 * (epsilon)/ewv) * (epsilon) - dmusig1 *
        sigmastar1)/ewv) * sigx1_1 * ewz * sigmastar1/(wzdeno *
      ewu1) + (((S^2 * ewu2 * (epsilon)^2/((sigma_sq2) *
      ewv) - 1) * suvq2 + ewu2/ssq2) * dmusig2 * ewu2/(sigma_sq2) +
      wusq2 * (S * ((2 * (S * dmusig2 * suvq2 * (epsilon)) +
        3 * pmusig2) * ewu2/(sigma_sq2) + S * wusq2 *
        sigx3_2 * (epsilon)/ewv) * (epsilon) - dmusig2 *
        sigmastar2)/ewv) * prC * sigx1_2 * sigmastar2/ewu2 -
      sigx5^2/sigx6)/sigx6, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx13_1 * (S * wusq1 *
      (epsilon)/ewv - sigx5/sigx6) + (((wusq1 * (pmusig1 -
      0.5 * (S * dmusig1 * ewu1 * (epsilon)/ssq1)) + S *
      sigx8_1 * ewu1 * (S * ewu1 * sigx2_1 * (epsilon)/(sigma_sq1) -
      2 * sigx3_1) * (epsilon)/sigmastar1) * sigmastar1 +
      0.5 * (wusq1 * ewu1 * ewv * sigx2_1/ssq1))/wzdeno -
      wzdeno * ewu1^2 * sigx2_1 * sigmastar1/(wzdeno *
        ewu1)^2)/(sigma_sq1)) * sigx1_1 * ewz/sigx6,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx13_2 * (S * wusq2 *
      (epsilon)/ewv - sigx5/sigx6) + ((wusq2 * (pmusig2 -
      0.5 * (S * dmusig2 * ewu2 * (epsilon)/ssq2)) + S *
      sigx8_2 * ewu2 * (S * ewu2 * sigx2_2 * (epsilon)/(sigma_sq2) -
      2 * sigx3_2) * (epsilon)/sigmastar2) * sigmastar2 +
      (0.5 * (wusq2 * ewu2 * ewv/ssq2) - sigmastar2) *
        sigx2_2)/(sigma_sq2)) * prC * sigx1_2/sigx6,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((sigx20_1 * sigx2_1 + (S * (0.5 * (wvsq1/ewv) +
    1/(sigma_sq1)) * dmusig1 * ewu1 * (epsilon)/sigmastar1 -
    pmusig1)/(sigma_sq1)) * ewu1/(sigma_sq1) + S * (2 * (sigx18_1 *
    ewu1^2/sigx19_1) - 4/(2 * ewv)^2) * sigx3_1 * (epsilon)) *
    sigmastar1 + 0.5 * (wvsq1 * ewu1^2 * sigx2_1/((sigma_sq1)^2 *
    sigmastar1)) + S * sigx24_1 * wusq1 * (epsilon)/ewv) *
    sigx1_1 * ewz/(wzdeno * ewu1) + (((sigx20_2 * sigx2_2 +
    (S * (0.5 * (wvsq2/ewv) + 1/(sigma_sq2)) * dmusig2 *
      ewu2 * (epsilon)/sigmastar2 - pmusig2)/(sigma_sq2)) *
    ewu2/(sigma_sq2) + S * (2 * (sigx18_2 * ewu2^2/sigx19_2) -
    4/(2 * ewv)^2) * sigx3_2 * (epsilon)) * sigmastar2 +
    0.5 * (wvsq2 * ewu2^2 * sigx2_2/((sigma_sq2)^2 * sigmastar2)) +
    S * sigx24_2 * wusq2 * (epsilon)/ewv) * prC * sigx1_2/ewu2 -
    sigx25 * sigx5/sigx6) * ewv/sigx6, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx23 * sigx1_1 * sigx4_1 *
      sigmastar1 - (sigx5 * sigx26/sigx6 + prC * sigx1_2 *
      sigx4_2 * sigmastar2/(wzdeno * ewu2))) * ewz/sigx6,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((wusq1 * (ewu1 * (S^2 * sigx8_1 * dmusig1 * (epsilon)^2 -
      sigx7_1/(sigma_sq1)) - 0.5 * (0.5 * (wusq1 * dmusig1 *
      ewv/sigmastar1) + S^2 * sigx8_1 * dmusig1 * ewu1 *
      (epsilon)^2)) + S^2 * ((sigx7_1 * wusq1 * sigx8_1/(sigma_sq1) -
      ((0.5 * (ewu1/(sigma_sq1)) + 1 - 0.5 * (0.5 * wusq1 +
        ewu1/(sigma_sq1))) * wusq1 * ewv/sigmastar1 +
        (2 - 2 * (s8sig1^2 * ewu1 * (sigma_sq1)/ssq1^2)) *
          sigmastar1) * sigx3_1/ssq1^2) * ewu1 + (1 -
      0.5 * wusq1) * sigx8_1 * sigx3_1) * ewu1 * (epsilon)^2/sigmastar1) *
      sigmastar1 + (0.5 * ((sigx7_1 * wusq1 + S * ewu1 *
      pmusig1 * (epsilon)/(sigma_sq1) - dmusig1 * sigmastar1) *
      ewu1/(sigma_sq1) - 0.5 * (wusq1 * sigx3_1)) + 0.5 *
      (sigx9_1 * ewu1/(sigma_sq1))) * wusq1 * ewv/sigmastar1)/sigx12_1 +
      ewu1 * (S^2 * sigx13_1 * sigx8_1 * ewu1 * (epsilon)^2/ssq1 -
        ((((sigx7_1 * wusq1 - S * pmusig1 * (epsilon)) *
          ewu1/(sigma_sq1) + dmusig1 * sigmastar1) *
          sigmastar1 + (0.5 * (wusq1 * ewv/ssq1) - 2 *
          (wzdeno^2 * ewu1 * sigmastar1/(wzdeno * ewu1)^2)) *
          sigx3_1 * ewu1)/(wzdeno * ewu1)^2 + sigx11_1/sigx12_1^2) *
          wzdeno) - sigx13_1^2 * sigx1_1 * ewz/sigx6) *
    sigx1_1 * ewz/sigx6, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx13_1 * sigx13_2 * prC * sigx1_1 *
      sigx1_2 * ewz/sigx6^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (wvsq1 * dmusig1) +
      0.5 * ((dmusig1 * ewv/(sigma_sq1) - S^2 * wvsq1 *
        sigx8_1 * dmusig1 * ewu1 * (epsilon)^2/sigmastar1) *
        ewu1/(sigma_sq1) - 0.5 * (wusq1 * wvsq1 * dmusig1)))/sigmastar1 +
      sigx7_1 * wusq1 * sigx20_1 + (S * (pmusig1 - ewu1 *
      (pmusig1/(sigma_sq1) + S * sigx8_1 * dmusig1 * (epsilon))) *
      (epsilon) - sigx15_1 * ewu1)/(sigma_sq1))/(sigma_sq1) -
      S^2 * (((0.5 * (wusq1 * ewv/(sigma_sq1)) + 0.5 *
        ((ewu1/(sigma_sq1) - 1) * ewv/(sigma_sq1) + 1 -
          0.5 * (wusq1 * wvsq1))) * ewu1/sigmastar1 +
        2 * sigx18_1)/sigx19_1 - ((ssq1^2 + 2 * (s8sig1 *
        (sigma_sq1)^2 * sigmastar1)) * sigmastar1 + 0.5 *
        (ssq1^2 * wusq1 * ewv/sigmastar1)) * sigx18_1 *
        ewu1/sigx19_1^2) * sigx3_1 * ewu1 * (epsilon)^2) *
      sigmastar1 + (0.5 * (sigx21_1 * wusq1 * ewv) + S^2 *
      sigx24_1 * sigx8_1 * ewu1 * (epsilon)^2)/ssq1 + 0.5 *
      (((sigx7_1 * wusq1 * wvsq1 + sigx3_1 * ewv/(sigma_sq1)) *
        ewu1/(sigma_sq1) + wvsq1 * sigx3_1)/ssq1 - s8sig1 *
        wvsq1 * sigx3_1 * ewu1/ssq1^2))/wzdeno - (sigx25 *
      sigx13_1/sigx6 + sigx24_1 * wzdeno * ewu1/(wzdeno *
      ewu1)^2)) * sigx1_1 * ewv * ewz/sigx6, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1 * sigx23/(sigma_sq1) -
      ((2 - 2 * (wzdeno^2 * ewu1^2/(wzdeno * ewu1)^2)) *
        ewz + 1) * sigx3_1/(wzdeno * ewu1)^2) * sigmastar1 +
      0.5 * (wusq1 * sigx23 * sigx3_1 * ewv/ssq1)) * ewu1 -
      sigx13_1 * sigx26 * ewz/sigx6) * sigx1_1 * ewz/sigx6,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((sigx7_2 * wusq2 +
      S * ewu2 * pmusig2 * (epsilon)/(sigma_sq2) - dmusig2 *
      sigmastar2) * ewu2/(sigma_sq2) - 0.5 * (wusq2 * sigx3_2)) +
      0.5 * (sigx9_2 * ewu2/(sigma_sq2)) - 0.5 * sigx3_2) *
      ewv/sigmastar2 - sigx7_2 * sigmastar2) * wusq2 +
      (wusq2 * (ewu2 * (S^2 * sigx8_2 * dmusig2 * (epsilon)^2 -
        sigx7_2/(sigma_sq2)) - 0.5 * (0.5 * (wusq2 *
        dmusig2 * ewv/sigmastar2) + S^2 * sigx8_2 * dmusig2 *
        ewu2 * (epsilon)^2)) + S^2 * ((sigx7_2 * wusq2 *
        sigx8_2/(sigma_sq2) - ((0.5 * (ewu2/(sigma_sq2)) +
        1 - 0.5 * (0.5 * wusq2 + ewu2/(sigma_sq2))) *
        wusq2 * ewv/sigmastar2 + (2 - 2 * (s8sig2^2 *
        ewu2 * (sigma_sq2)/ssq2^2)) * sigmastar2) * sigx3_2/ssq2^2) *
        ewu2 + (1 - 0.5 * wusq2) * sigx8_2 * sigx3_2) *
        ewu2 * (epsilon)^2/sigmastar2) * sigmastar2 +
      ewu2 * (S^2 * sigx13_2 * sigx8_2 * ewu2 * (epsilon)^2/sigmastar2 -
        sigx11_2/(sigma_sq2)))/(sigma_sq2) + sigx3_2 *
      sigmastar2/ewu2 - sigx13_2^2 * prC * sigx1_2/sigx6) *
      prC * sigx1_2/sigx6, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (wvsq2 * dmusig2) +
      0.5 * ((dmusig2 * ewv/(sigma_sq2) - S^2 * wvsq2 *
        sigx8_2 * dmusig2 * ewu2 * (epsilon)^2/sigmastar2) *
        ewu2/(sigma_sq2) - 0.5 * (wusq2 * wvsq2 * dmusig2)))/sigmastar2 +
      sigx7_2 * wusq2 * sigx20_2 + (S * (pmusig2 - ewu2 *
      (pmusig2/(sigma_sq2) + S * sigx8_2 * dmusig2 * (epsilon))) *
      (epsilon) - sigx15_2 * ewu2)/(sigma_sq2))/(sigma_sq2) -
      S^2 * (((0.5 * (wusq2 * ewv/(sigma_sq2)) + 0.5 *
        ((ewu2/(sigma_sq2) - 1) * ewv/(sigma_sq2) + 1 -
          0.5 * (wusq2 * wvsq2))) * ewu2/sigmastar2 +
        2 * sigx18_2)/sigx19_2 - ((ssq2^2 + 2 * (s8sig2 *
        (sigma_sq2)^2 * sigmastar2)) * sigmastar2 + 0.5 *
        (ssq2^2 * wusq2 * ewv/sigmastar2)) * sigx18_2 *
        ewu2/sigx19_2^2) * sigx3_2 * ewu2 * (epsilon)^2) *
      sigmastar2 + (0.5 * (sigx21_2 * wusq2 * ewv) + S^2 *
      sigx24_2 * sigx8_2 * ewu2 * (epsilon)^2)/ssq2 + 0.5 *
      (((sigx7_2 * wusq2 * wvsq2 + sigx3_2 * ewv/(sigma_sq2)) *
        ewu2/(sigma_sq2) + wvsq2 * sigx3_2)/ssq2 - s8sig2 *
        wvsq2 * sigx3_2 * ewu2/ssq2^2) - (sigx25 * sigx13_2/sigx6 +
      sigx24_2/ewu2)) * prC * sigx1_2 * ewv/sigx6, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx13_2 * sigx26/sigx6 + sigx11_2/(wzdeno * (sigma_sq2)) -
      wzdeno * sigx3_2 * ewu2 * sigmastar2/(wzdeno * ewu2)^2) *
      prC * sigx1_2 * ewz/sigx6), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((sigx15_1 * sigx20_1 +
      (S * (S * sigx18_1 * dmusig1 * ewu1 * (epsilon)/ssq1^2 -
        2 * (pmusig1/(sigma_sq1))) * (epsilon) - 0.5 *
        sigx14_1)/(sigma_sq1)) * ewv + (0.5 * (ewv *
      (S^2 * sigx18_1 * dmusig1 * ewu1^2 * (epsilon)^2/(ssq1^2 *
        sigmastar1) - dmusig1)/(sigma_sq1) - 0.5 * (wvsq1 *
      dmusig1)) + 0.5 * dmusig1) * wvsq1/sigmastar1 + S *
      pmusig1 * (epsilon)/(sigma_sq1)) * ewu1/(sigma_sq1) +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx18_1 * (1/sigx19_1 - ((ssq1^2 +
        2 * (sigx18_1 * (sigma_sq1)^2 * sigmastar1)) *
        sigmastar1 + 0.5 * (ssq1^2 * wvsq1 * ewu1/sigmastar1)) *
        ewv/sigx19_1^2) + (0.5 * (ewv/(sigma_sq1)) -
        0.5 * (0.5 * wvsq1 + ewv/(sigma_sq1))) * wvsq1/(ssq1^2 *
        ewv)) * ewu1^2 * (epsilon)^2) * sigx3_1) * sigmastar1 +
      (sigx24_1 * sigx20_1 + (0.5 * (((0.5 * sigx14_1 +
        2 * (S * pmusig1 * (epsilon)/(sigma_sq1))) *
        ewu1 - dmusig1 * sigmastar1)/((sigma_sq1)^2 *
        sigmastar1) - sigx18_1 * sigx3_1/ssq1^2) + 0.5 *
        (sigx21_1/ssq1)) * wvsq1 * ewu1) * ewv + 0.5 *
      sigx22_1) * sigx1_1 * ewz/(wzdeno * ewu1) + ((((sigx15_2 *
      sigx20_2 + (S * (S * sigx18_2 * dmusig2 * ewu2 *
      (epsilon)/ssq2^2 - 2 * (pmusig2/(sigma_sq2))) * (epsilon) -
      0.5 * sigx14_2)/(sigma_sq2)) * ewv + (0.5 * (ewv *
      (S^2 * sigx18_2 * dmusig2 * ewu2^2 * (epsilon)^2/(ssq2^2 *
        sigmastar2) - dmusig2)/(sigma_sq2) - 0.5 * (wvsq2 *
      dmusig2)) + 0.5 * dmusig2) * wvsq2/sigmastar2 + S *
      pmusig2 * (epsilon)/(sigma_sq2)) * ewu2/(sigma_sq2) +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx18_2 * (1/sigx19_2 - ((ssq2^2 +
        2 * (sigx18_2 * (sigma_sq2)^2 * sigmastar2)) *
        sigmastar2 + 0.5 * (ssq2^2 * wvsq2 * ewu2/sigmastar2)) *
        ewv/sigx19_2^2) + (0.5 * (ewv/(sigma_sq2)) -
        0.5 * (0.5 * wvsq2 + ewv/(sigma_sq2))) * wvsq2/(ssq2^2 *
        ewv)) * ewu2^2 * (epsilon)^2) * sigx3_2) * sigmastar2 +
      (sigx24_2 * sigx20_2 + (0.5 * (((0.5 * sigx14_2 +
        2 * (S * pmusig2 * (epsilon)/(sigma_sq2))) *
        ewu2 - dmusig2 * sigmastar2)/((sigma_sq2)^2 *
        sigmastar2) - sigx18_2 * sigx3_2/ssq2^2) + 0.5 *
        (sigx21_2/ssq2)) * wvsq2 * ewu2) * ewv + 0.5 *
      sigx22_2) * prC * sigx1_2/ewu2 - sigx25^2 * ewv/sigx6) *
      ewv/sigx6, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (sigx24_1 * sigx23 * sigx1_1 - (sigx25 *
      sigx26/sigx6 + sigx24_2 * prC * sigx1_2/(wzdeno *
      ewu2))) * ewv * ewz/sigx6, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewu2) +
      ewu2/(wzdeno * ewu2)^2) * sigx3_2 * sigx1_2 * sigmastar2 -
      (sigx26^2/sigx6 + (2 - 2 * (wzdeno * ewu1^2 * ewz/(wzdeno *
        ewu1)^2)) * sigx3_1 * sigx1_1 * ewu1 * sigmastar1/(wzdeno *
        ewu1)^2)) * ewz + sigx23 * sigx3_1 * sigx1_1 *
      sigmastar1 - prC * sigx3_2 * sigx1_2 * sigmastar2/(wzdeno *
      ewu2)) * ewz/sigx6, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cauchit specification class membership
chessmisfraynormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  musig1 <- (S * ewu1 * (epsilon)/starsq1)
  musig2 <- (S * ewu2 * (epsilon)/starsq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  wuv1 <- (sigmastar1/ewv - ewu1/starsq1)
  wuv2 <- (sigmastar2/ewv - ewu2/starsq2)
  dvs1 <- (dmusig1 * ewv/sigmastar1)
  dvs2 <- (dmusig2 * ewv/sigmastar2)
  dpsv1 <- (0.5 * dvs1 - S * pmusig1 * (epsilon))
  dpsv2 <- (0.5 * dvs2 - S * pmusig2 * (epsilon))
  vepsi <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx2_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/sigma_sq1)
  sigx2_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/sigma_sq2)
  sigx3_1 <- (pmusig1 + S * dmusig1 * wuv1 * (epsilon))
  sigx3_2 <- (pmusig2 + S * dmusig2 * wuv2 * (epsilon))
  sigx4_1 <- (ewu1 * sigx3_1/sigma_sq1 + S * usq1 * sigx2_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx3_2/sigma_sq2 + S * usq2 * sigx2_2 *
    (epsilon)/ewv)
  sigx5_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx5_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx6_1 <- (1/starsq1 - sigx5_1 * ewu1/starsq1^2)
  sigx6_2 <- (1/starsq2 - sigx5_2 * ewu2/starsq2^2)
  sigx7_1 <- (dpsv1 * usq1 + S^2 * sigx6_1 * sigx2_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx7_2 <- (dpsv2 * usq2 + S^2 * sigx6_2 * sigx2_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx8_1 <- (usq1 * sigx2_1 * ewv/sigmastar1)
  sigx8_2 <- (usq2 * sigx2_2 * ewv/sigmastar2)
  sigx9_1 <- (sigx7_1 * sigmastar1 + 0.5 * sigx8_1)
  sigx9_2 <- (sigx7_2 * sigmastar2 + 0.5 * sigx8_2)
  sigx10_1 <- (vsq1 * dmusig1/sigmastar1)
  sigx10_2 <- (vsq2 * dmusig2/sigmastar2)
  sigx11_1 <- (0.5 * sigx10_1 + S * pmusig1 * (epsilon)/sigma_sq1)
  sigx11_2 <- (0.5 * sigx10_2 + S * pmusig2 * (epsilon)/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (starsq1^2 * sigma_sq1 * sigmastar1)
  sigx13_2 <- (starsq2^2 * sigma_sq2 * sigmastar2)
  sigx14_1 <- (2 * vepsi - S^2 * sigx12_1 * ewu1^2 * (epsilon)^2/sigx13_1)
  sigx14_2 <- (2 * vepsi - S^2 * sigx12_2 * ewu2^2 * (epsilon)^2/sigx13_2)
  sigx15_1 <- (sigx11_1 * ewu1/sigma_sq1 + sigx14_1 * sigx2_1)
  sigx15_2 <- (sigx11_2 * ewu2/sigma_sq2 + sigx14_2 * sigx2_2)
  sigx16_1 <- (sigx15_1 * sigmastar1 + 0.5 * (vsq1 * sigx2_1 *
    ewu1/starsq1))
  sigx16_2 <- (sigx15_2 * sigmastar2 + 0.5 * (vsq2 * sigx2_2 *
    ewu2/starsq2))
  sigx17_1 <- (sigx9_1/sigma_sq1 - sigx2_1 * sigmastar1/ewu1)
  sigx17_2 <- (sigx9_2/sigma_sq2 - sigx2_2 * sigmastar2/ewu2)
  sigx18 <- (ewz2 * sigx1_2 * sigx4_2 * sigmastar2/ewu2 + ewz1 *
    sigx1_1 * sigx4_1 * sigmastar1/ewu1)
  sigx19 <- (ewz2 * sigx2_2 * sigx1_2 * sigmastar2/ewu2 + ewz1 *
    sigx2_1 * sigx1_1 * sigmastar1/ewu1)
  sigx20 <- (sigx16_1 * ewz1 * sigx1_1/ewu1 + sigx16_2 * ewz2 *
    sigx1_2/ewu2)
  sigx21 <- (sigx2_1 * sigx1_1 * sigmastar1/ewu1 - sigx2_2 *
    sigx1_2 * sigmastar2/ewu2)
  sigx22 <- (pi * sigx19 * ((Wz)^2 + 1))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((S^2 * ewu1 * (epsilon)^2/(sigma_sq1 *
      ewv) - 1) * wuv1 + ewu1/starsq1) * dmusig1 * ewu1/sigma_sq1 +
      usq1 * (S * ((2 * (S * dmusig1 * wuv1 * (epsilon)) +
        3 * pmusig1) * ewu1/sigma_sq1 + S * usq1 * sigx2_1 *
        (epsilon)/ewv) * (epsilon) - dmusig1 * sigmastar1)/ewv) *
      ewz1 * sigx1_1 * sigmastar1/ewu1 + (((S^2 * ewu2 *
      (epsilon)^2/(sigma_sq2 * ewv) - 1) * wuv2 + ewu2/starsq2) *
      dmusig2 * ewu2/sigma_sq2 + usq2 * (S * ((2 * (S *
      dmusig2 * wuv2 * (epsilon)) + 3 * pmusig2) * ewu2/sigma_sq2 +
      S * usq2 * sigx2_2 * (epsilon)/ewv) * (epsilon) -
      dmusig2 * sigmastar2)/ewv) * ewz2 * sigx1_2 * sigmastar2/ewu2 -
      sigx18^2/sigx19)/sigx19, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx17_1 * (S * usq1 *
      (epsilon)/ewv - sigx18/sigx19) + ((usq1 * (pmusig1 -
      0.5 * (S * dmusig1 * ewu1 * (epsilon)/starsq1)) +
      S * sigx6_1 * ewu1 * (S * ewu1 * sigx3_1 * (epsilon)/sigma_sq1 -
        2 * sigx2_1) * (epsilon)/sigmastar1) * sigmastar1 +
      (0.5 * (usq1 * ewu1 * ewv/starsq1) - sigmastar1) *
        sigx3_1)/sigma_sq1) * ewz1 * sigx1_1/sigx19,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx17_2 * (S * usq2 *
      (epsilon)/ewv - sigx18/sigx19) + ((usq2 * (pmusig2 -
      0.5 * (S * dmusig2 * ewu2 * (epsilon)/starsq2)) +
      S * sigx6_2 * ewu2 * (S * ewu2 * sigx3_2 * (epsilon)/sigma_sq2 -
        2 * sigx2_2) * (epsilon)/sigmastar2) * sigmastar2 +
      (0.5 * (usq2 * ewu2 * ewv/starsq2) - sigmastar2) *
        sigx3_2)/sigma_sq2) * ewz2 * sigx1_2/sigx19,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((sigx14_1 * sigx3_1 + (S * (0.5 * (vsq1/ewv) +
    1/sigma_sq1) * dmusig1 * ewu1 * (epsilon)/sigmastar1 -
    pmusig1)/sigma_sq1) * ewu1/sigma_sq1 + S * (2 * (sigx12_1 *
    ewu1^2/sigx13_1) - 4/(2 * ewv)^2) * sigx2_1 * (epsilon)) *
    sigmastar1 + 0.5 * (vsq1 * ewu1^2 * sigx3_1/(sigma_sq1^2 *
    sigmastar1)) + S * sigx16_1 * usq1 * (epsilon)/ewv) *
    ewz1 * sigx1_1/ewu1 + (((sigx14_2 * sigx3_2 + (S * (0.5 *
    (vsq2/ewv) + 1/sigma_sq2) * dmusig2 * ewu2 * (epsilon)/sigmastar2 -
    pmusig2)/sigma_sq2) * ewu2/sigma_sq2 + S * (2 * (sigx12_2 *
    ewu2^2/sigx13_2) - 4/(2 * ewv)^2) * sigx2_2 * (epsilon)) *
    sigmastar2 + 0.5 * (vsq2 * ewu2^2 * sigx3_2/(sigma_sq2^2 *
    sigmastar2)) + S * sigx16_2 * usq2 * (epsilon)/ewv) *
    ewz2 * sigx1_2/ewu2 - sigx20 * sigx18/sigx19) * ewv/sigx19,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx1_1 * sigx4_1 *
      sigmastar1/ewu1 - sigx1_2 * sigx4_2 * sigmastar2/ewu2)/sigx22 -
      pi * sigx18 * sigx21 * ((Wz)^2 + 1)/sigx22^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((dpsv1 * usq1 + S * ewu1 * pmusig1 * (epsilon)/sigma_sq1 -
      dmusig1 * sigmastar1) * ewu1/sigma_sq1 - 0.5 * (usq1 *
      sigx2_1)) + 0.5 * (sigx7_1 * ewu1/sigma_sq1) - 0.5 *
      sigx2_1) * ewv/sigmastar1 - dpsv1 * sigmastar1) *
      usq1 + (usq1 * (ewu1 * (S^2 * sigx6_1 * dmusig1 *
      (epsilon)^2 - dpsv1/sigma_sq1) - 0.5 * (0.5 * (usq1 *
      dmusig1 * ewv/sigmastar1) + S^2 * sigx6_1 * dmusig1 *
      ewu1 * (epsilon)^2)) + S^2 * ((dpsv1 * usq1 * sigx6_1/sigma_sq1 -
      ((0.5 * (ewu1/sigma_sq1) + 1 - 0.5 * (0.5 * usq1 +
        ewu1/sigma_sq1)) * usq1 * ewv/sigmastar1 + (2 -
        2 * (sigx5_1^2 * ewu1 * sigma_sq1/starsq1^2)) *
        sigmastar1) * sigx2_1/starsq1^2) * ewu1 + (1 -
      0.5 * usq1) * sigx6_1 * sigx2_1) * ewu1 * (epsilon)^2/sigmastar1) *
      sigmastar1 + ewu1 * (S^2 * sigx17_1 * sigx6_1 * ewu1 *
      (epsilon)^2/sigmastar1 - sigx9_1/sigma_sq1))/sigma_sq1 +
      sigx2_1 * sigmastar1/ewu1 - sigx17_1^2 * ewz1 * sigx1_1/sigx19) *
    ewz1 * sigx1_1/sigx19, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx17_1 * sigx17_2 * ewz2 * ewz1 *
      sigx1_1 * sigx1_2/sigx19^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (vsq1 * dmusig1) +
      0.5 * ((dmusig1 * ewv/sigma_sq1 - S^2 * vsq1 * sigx6_1 *
        dmusig1 * ewu1 * (epsilon)^2/sigmastar1) * ewu1/sigma_sq1 -
        0.5 * (usq1 * vsq1 * dmusig1)))/sigmastar1 +
      dpsv1 * usq1 * sigx14_1 + (S * (pmusig1 - ewu1 *
      (pmusig1/sigma_sq1 + S * sigx6_1 * dmusig1 * (epsilon))) *
      (epsilon) - sigx11_1 * ewu1)/sigma_sq1)/sigma_sq1 -
      S^2 * (((0.5 * (usq1 * ewv/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
        1) * ewv/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) *
        ewu1/sigmastar1 + 2 * sigx12_1)/sigx13_1 - ((starsq1^2 +
        2 * (sigx5_1 * sigma_sq1^2 * sigmastar1)) * sigmastar1 +
        0.5 * (starsq1^2 * usq1 * ewv/sigmastar1)) *
        sigx12_1 * ewu1/sigx13_1^2) * sigx2_1 * ewu1 *
        (epsilon)^2) * sigmastar1 + (0.5 * (sigx15_1 *
      usq1 * ewv) + S^2 * sigx16_1 * sigx6_1 * ewu1 * (epsilon)^2)/starsq1 +
      0.5 * (((dpsv1 * usq1 * vsq1 + sigx2_1 * ewv/sigma_sq1) *
        ewu1/sigma_sq1 + vsq1 * sigx2_1)/starsq1 - sigx5_1 *
        vsq1 * sigx2_1 * ewu1/starsq1^2) - (sigx20 *
      sigx17_1/sigx19 + sigx16_1/ewu1)) * ewz1 * sigx1_1 *
      ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx17_1 * (1/sigx22 - pi *
      sigx21 * ((Wz)^2 + 1) * ewz1/sigx22^2) * sigx1_1,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((dpsv2 * usq2 +
      S * ewu2 * pmusig2 * (epsilon)/sigma_sq2 - dmusig2 *
      sigmastar2) * ewu2/sigma_sq2 - 0.5 * (usq2 * sigx2_2)) +
      0.5 * (sigx7_2 * ewu2/sigma_sq2) - 0.5 * sigx2_2) *
      ewv/sigmastar2 - dpsv2 * sigmastar2) * usq2 + (usq2 *
      (ewu2 * (S^2 * sigx6_2 * dmusig2 * (epsilon)^2 -
        dpsv2/sigma_sq2) - 0.5 * (0.5 * (usq2 * dmusig2 *
        ewv/sigmastar2) + S^2 * sigx6_2 * dmusig2 * ewu2 *
        (epsilon)^2)) + S^2 * ((dpsv2 * usq2 * sigx6_2/sigma_sq2 -
      ((0.5 * (ewu2/sigma_sq2) + 1 - 0.5 * (0.5 * usq2 +
        ewu2/sigma_sq2)) * usq2 * ewv/sigmastar2 + (2 -
        2 * (sigx5_2^2 * ewu2 * sigma_sq2/starsq2^2)) *
        sigmastar2) * sigx2_2/starsq2^2) * ewu2 + (1 -
      0.5 * usq2) * sigx6_2 * sigx2_2) * ewu2 * (epsilon)^2/sigmastar2) *
      sigmastar2 + ewu2 * (S^2 * sigx17_2 * sigx6_2 * ewu2 *
      (epsilon)^2/sigmastar2 - sigx9_2/sigma_sq2))/sigma_sq2 +
      sigx2_2 * sigmastar2/ewu2 - sigx17_2^2 * ewz2 * sigx1_2/sigx19) *
      ewz2 * sigx1_2/sigx19, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (vsq2 * dmusig2) +
      0.5 * ((dmusig2 * ewv/sigma_sq2 - S^2 * vsq2 * sigx6_2 *
        dmusig2 * ewu2 * (epsilon)^2/sigmastar2) * ewu2/sigma_sq2 -
        0.5 * (usq2 * vsq2 * dmusig2)))/sigmastar2 +
      dpsv2 * usq2 * sigx14_2 + (S * (pmusig2 - ewu2 *
      (pmusig2/sigma_sq2 + S * sigx6_2 * dmusig2 * (epsilon))) *
      (epsilon) - sigx11_2 * ewu2)/sigma_sq2)/sigma_sq2 -
      S^2 * (((0.5 * (usq2 * ewv/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
        1) * ewv/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
        ewu2/sigmastar2 + 2 * sigx12_2)/sigx13_2 - ((starsq2^2 +
        2 * (sigx5_2 * sigma_sq2^2 * sigmastar2)) * sigmastar2 +
        0.5 * (starsq2^2 * usq2 * ewv/sigmastar2)) *
        sigx12_2 * ewu2/sigx13_2^2) * sigx2_2 * ewu2 *
        (epsilon)^2) * sigmastar2 + (0.5 * (sigx15_2 *
      usq2 * ewv) + S^2 * sigx16_2 * sigx6_2 * ewu2 * (epsilon)^2)/starsq2 +
      0.5 * (((dpsv2 * usq2 * vsq2 + sigx2_2 * ewv/sigma_sq2) *
        ewu2/sigma_sq2 + vsq2 * sigx2_2)/starsq2 - sigx5_2 *
        vsq2 * sigx2_2 * ewu2/starsq2^2) - (sigx20 *
      sigx17_2/sigx19 + sigx16_2/ewu2)) * ewz2 * sigx1_2 *
      ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx17_2 * (1/sigx22 + pi * sigx21 * ((Wz)^2 + 1) *
      ewz2/sigx22^2) * sigx1_2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((sigx11_1 * sigx14_1 +
      (S * (S * sigx12_1 * dmusig1 * ewu1 * (epsilon)/starsq1^2 -
        2 * (pmusig1/sigma_sq1)) * (epsilon) - 0.5 *
        sigx10_1)/sigma_sq1) * ewv + (0.5 * (ewv * (S^2 *
      sigx12_1 * dmusig1 * ewu1^2 * (epsilon)^2/(starsq1^2 *
      sigmastar1) - dmusig1)/sigma_sq1 - 0.5 * (vsq1 *
      dmusig1)) + 0.5 * dmusig1) * vsq1/sigmastar1 + S *
      pmusig1 * (epsilon)/sigma_sq1) * ewu1/sigma_sq1 +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx12_1 * (1/sigx13_1 - ((starsq1^2 +
        2 * (sigx12_1 * sigma_sq1^2 * sigmastar1)) *
        sigmastar1 + 0.5 * (starsq1^2 * vsq1 * ewu1/sigmastar1)) *
        ewv/sigx13_1^2) + (0.5 * (ewv/sigma_sq1) - 0.5 *
        (0.5 * vsq1 + ewv/sigma_sq1)) * vsq1/(starsq1^2 *
        ewv)) * ewu1^2 * (epsilon)^2) * sigx2_1) * sigmastar1 +
      (sigx16_1 * sigx14_1 + (0.5 * (((0.5 * sigx10_1 +
        2 * (S * pmusig1 * (epsilon)/sigma_sq1)) * ewu1 -
        dmusig1 * sigmastar1)/(sigma_sq1^2 * sigmastar1) -
        sigx12_1 * sigx2_1/starsq1^2) + 0.5 * (sigx15_1/starsq1)) *
        vsq1 * ewu1) * ewv + 0.5 * (vsq1 * sigx2_1 *
      ewu1/starsq1)) * ewz1 * sigx1_1/ewu1 + ((((sigx11_2 *
      sigx14_2 + (S * (S * sigx12_2 * dmusig2 * ewu2 *
      (epsilon)/starsq2^2 - 2 * (pmusig2/sigma_sq2)) *
      (epsilon) - 0.5 * sigx10_2)/sigma_sq2) * ewv + (0.5 *
      (ewv * (S^2 * sigx12_2 * dmusig2 * ewu2^2 * (epsilon)^2/(starsq2^2 *
        sigmastar2) - dmusig2)/sigma_sq2 - 0.5 * (vsq2 *
        dmusig2)) + 0.5 * dmusig2) * vsq2/sigmastar2 +
      S * pmusig2 * (epsilon)/sigma_sq2) * ewu2/sigma_sq2 +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx12_2 * (1/sigx13_2 - ((starsq2^2 +
        2 * (sigx12_2 * sigma_sq2^2 * sigmastar2)) *
        sigmastar2 + 0.5 * (starsq2^2 * vsq2 * ewu2/sigmastar2)) *
        ewv/sigx13_2^2) + (0.5 * (ewv/sigma_sq2) - 0.5 *
        (0.5 * vsq2 + ewv/sigma_sq2)) * vsq2/(starsq2^2 *
        ewv)) * ewu2^2 * (epsilon)^2) * sigx2_2) * sigmastar2 +
      (sigx16_2 * sigx14_2 + (0.5 * (((0.5 * sigx10_2 +
        2 * (S * pmusig2 * (epsilon)/sigma_sq2)) * ewu2 -
        dmusig2 * sigmastar2)/(sigma_sq2^2 * sigmastar2) -
        sigx12_2 * sigx2_2/starsq2^2) + 0.5 * (sigx15_2/starsq2)) *
        vsq2 * ewu2) * ewv + 0.5 * (vsq2 * sigx2_2 *
      ewu2/starsq2)) * ewz2 * sigx1_2/ewu2 - sigx20^2 *
      ewv/sigx19) * ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx16_1 * sigx1_1/ewu1 - sigx16_2 *
      sigx1_2/ewu2)/sigx22 - pi * sigx20 * sigx21 * ((Wz)^2 +
      1)/sigx22^2) * ewv, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx21 * (sigx2_1 * sigx1_1 *
      sigmastar1/ewu1 + 2 * (pi * Wz * sigx19) - sigx2_2 *
      sigx1_2 * sigmastar2/ewu2)/sigx22^2), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# probit specification class membership
chessmisfraynormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  musig1 <- (S * ewu1 * (epsilon)/starsq1)
  musig2 <- (S * ewu2 * (epsilon)/starsq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  wuv1 <- (sigmastar1/ewv - ewu1/starsq1)
  wuv2 <- (sigmastar2/ewv - ewu2/starsq2)
  dvs1 <- (dmusig1 * ewv/sigmastar1)
  dvs2 <- (dmusig2 * ewv/sigmastar2)
  dpsv1 <- (0.5 * dvs1 - S * pmusig1 * (epsilon))
  dpsv2 <- (0.5 * dvs2 - S * pmusig2 * (epsilon))
  vepsi <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx2_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/sigma_sq1)
  sigx2_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/sigma_sq2)
  sigx3_1 <- (pmusig1 + S * dmusig1 * wuv1 * (epsilon))
  sigx3_2 <- (pmusig2 + S * dmusig2 * wuv2 * (epsilon))
  sigx4_1 <- (ewu1 * sigx3_1/sigma_sq1 + S * usq1 * sigx2_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx3_2/sigma_sq2 + S * usq2 * sigx2_2 *
    (epsilon)/ewv)
  sigx5_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx5_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx6_1 <- (1/starsq1 - sigx5_1 * ewu1/starsq1^2)
  sigx6_2 <- (1/starsq2 - sigx5_2 * ewu2/starsq2^2)
  sigx7_1 <- (dpsv1 * usq1 + S^2 * sigx6_1 * sigx2_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx7_2 <- (dpsv2 * usq2 + S^2 * sigx6_2 * sigx2_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx8_1 <- (usq1 * sigx2_1 * ewv/sigmastar1)
  sigx8_2 <- (usq2 * sigx2_2 * ewv/sigmastar2)
  sigx9_1 <- (sigx7_1 * sigmastar1 + 0.5 * sigx8_1)
  sigx9_2 <- (sigx7_2 * sigmastar2 + 0.5 * sigx8_2)
  sigx10_1 <- (vsq1 * dmusig1/sigmastar1)
  sigx10_2 <- (vsq2 * dmusig2/sigmastar2)
  sigx11_1 <- (0.5 * sigx10_1 + S * pmusig1 * (epsilon)/sigma_sq1)
  sigx11_2 <- (0.5 * sigx10_2 + S * pmusig2 * (epsilon)/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (starsq1^2 * sigma_sq1 * sigmastar1)
  sigx13_2 <- (starsq2^2 * sigma_sq2 * sigmastar2)
  sigx14_1 <- (2 * vepsi - S^2 * sigx12_1 * ewu1^2 * (epsilon)^2/sigx13_1)
  sigx14_2 <- (2 * vepsi - S^2 * sigx12_2 * ewu2^2 * (epsilon)^2/sigx13_2)
  sigx15_1 <- (sigx11_1 * ewu1/sigma_sq1 + sigx14_1 * sigx2_1)
  sigx15_2 <- (sigx11_2 * ewu2/sigma_sq2 + sigx14_2 * sigx2_2)
  sigx16_1 <- (sigx15_1 * sigmastar1 + 0.5 * (vsq1 * sigx2_1 *
    ewu1/starsq1))
  sigx16_2 <- (sigx15_2 * sigmastar2 + 0.5 * (vsq2 * sigx2_2 *
    ewu2/starsq2))
  sigx17_1 <- (sigx9_1/sigma_sq1 - sigx2_1 * sigmastar1/ewu1)
  sigx17_2 <- (sigx9_2/sigma_sq2 - sigx2_2 * sigmastar2/ewu2)
  sigx18 <- ((1 - pwZ) * sigx1_2 * sigx4_2 * sigmastar2/ewu2 +
    sigx1_1 * sigx4_1 * pwZ * sigmastar1/ewu1)
  sigx19 <- ((1 - pwZ) * sigx2_2 * sigx1_2 * sigmastar2/ewu2 +
    sigx2_1 * sigx1_1 * pwZ * sigmastar1/ewu1)
  sigx20 <- (sigx16_1 * sigx1_1 * pwZ/ewu1 + sigx16_2 * (1 -
    pwZ) * sigx1_2/ewu2)
  sigx21 <- (sigx2_1 * sigx1_1 * sigmastar1/ewu1 - sigx2_2 *
    sigx1_2 * sigmastar2/ewu2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((S^2 * ewu1 * (epsilon)^2/(sigma_sq1 *
      ewv) - 1) * wuv1 + ewu1/starsq1) * dmusig1 * ewu1/sigma_sq1 +
      usq1 * (S * ((2 * (S * dmusig1 * wuv1 * (epsilon)) +
        3 * pmusig1) * ewu1/sigma_sq1 + S * usq1 * sigx2_1 *
        (epsilon)/ewv) * (epsilon) - dmusig1 * sigmastar1)/ewv) *
      pwZ * sigx1_1 * sigmastar1/ewu1 + (((S^2 * ewu2 *
      (epsilon)^2/(sigma_sq2 * ewv) - 1) * wuv2 + ewu2/starsq2) *
      dmusig2 * ewu2/sigma_sq2 + usq2 * (S * ((2 * (S *
      dmusig2 * wuv2 * (epsilon)) + 3 * pmusig2) * ewu2/sigma_sq2 +
      S * usq2 * sigx2_2 * (epsilon)/ewv) * (epsilon) -
      dmusig2 * sigmastar2)/ewv) * (1 - pwZ) * sigx1_2 *
      sigmastar2/ewu2 - sigx18^2/sigx19)/sigx19, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx17_1 * (S * usq1 *
      (epsilon)/ewv - sigx18/sigx19) + ((usq1 * (pmusig1 -
      0.5 * (S * dmusig1 * ewu1 * (epsilon)/starsq1)) +
      S * sigx6_1 * ewu1 * (S * ewu1 * sigx3_1 * (epsilon)/sigma_sq1 -
        2 * sigx2_1) * (epsilon)/sigmastar1) * sigmastar1 +
      (0.5 * (usq1 * ewu1 * ewv/starsq1) - sigmastar1) *
        sigx3_1)/sigma_sq1) * pwZ * sigx1_1/sigx19, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx17_2 * (S * usq2 *
      (epsilon)/ewv - sigx18/sigx19) + ((usq2 * (pmusig2 -
      0.5 * (S * dmusig2 * ewu2 * (epsilon)/starsq2)) +
      S * sigx6_2 * ewu2 * (S * ewu2 * sigx3_2 * (epsilon)/sigma_sq2 -
        2 * sigx2_2) * (epsilon)/sigmastar2) * sigmastar2 +
      (0.5 * (usq2 * ewu2 * ewv/starsq2) - sigmastar2) *
        sigx3_2)/sigma_sq2) * (1 - pwZ) * sigx1_2/sigx19,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((sigx14_1 * sigx3_1 + (S * (0.5 * (vsq1/ewv) +
    1/sigma_sq1) * dmusig1 * ewu1 * (epsilon)/sigmastar1 -
    pmusig1)/sigma_sq1) * ewu1/sigma_sq1 + S * (2 * (sigx12_1 *
    ewu1^2/sigx13_1) - 4/(2 * ewv)^2) * sigx2_1 * (epsilon)) *
    sigmastar1 + 0.5 * (vsq1 * ewu1^2 * sigx3_1/(sigma_sq1^2 *
    sigmastar1)) + S * sigx16_1 * usq1 * (epsilon)/ewv) *
    pwZ * sigx1_1/ewu1 + (((sigx14_2 * sigx3_2 + (S * (0.5 *
    (vsq2/ewv) + 1/sigma_sq2) * dmusig2 * ewu2 * (epsilon)/sigmastar2 -
    pmusig2)/sigma_sq2) * ewu2/sigma_sq2 + S * (2 * (sigx12_2 *
    ewu2^2/sigx13_2) - 4/(2 * ewv)^2) * sigx2_2 * (epsilon)) *
    sigmastar2 + 0.5 * (vsq2 * ewu2^2 * sigx3_2/(sigma_sq2^2 *
    sigmastar2)) + S * sigx16_2 * usq2 * (epsilon)/ewv) *
    (1 - pwZ) * sigx1_2/ewu2 - sigx20 * sigx18/sigx19) *
    ewv/sigx19, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * dwZ * (sigx1_1 * sigx4_1 *
      sigmastar1/ewu1 - (sigx18 * sigx21/sigx19 + sigx1_2 *
      sigx4_2 * sigmastar2/ewu2))/sigx19, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((dpsv1 * usq1 + S * ewu1 * pmusig1 * (epsilon)/sigma_sq1 -
      dmusig1 * sigmastar1) * ewu1/sigma_sq1 - 0.5 * (usq1 *
      sigx2_1)) + 0.5 * (sigx7_1 * ewu1/sigma_sq1) - 0.5 *
      sigx2_1) * ewv/sigmastar1 - dpsv1 * sigmastar1) *
      usq1 + (usq1 * (ewu1 * (S^2 * sigx6_1 * dmusig1 *
      (epsilon)^2 - dpsv1/sigma_sq1) - 0.5 * (0.5 * (usq1 *
      dmusig1 * ewv/sigmastar1) + S^2 * sigx6_1 * dmusig1 *
      ewu1 * (epsilon)^2)) + S^2 * ((dpsv1 * usq1 * sigx6_1/sigma_sq1 -
      ((0.5 * (ewu1/sigma_sq1) + 1 - 0.5 * (0.5 * usq1 +
        ewu1/sigma_sq1)) * usq1 * ewv/sigmastar1 + (2 -
        2 * (sigx5_1^2 * ewu1 * sigma_sq1/starsq1^2)) *
        sigmastar1) * sigx2_1/starsq1^2) * ewu1 + (1 -
      0.5 * usq1) * sigx6_1 * sigx2_1) * ewu1 * (epsilon)^2/sigmastar1) *
      sigmastar1 + ewu1 * (S^2 * sigx17_1 * sigx6_1 * ewu1 *
      (epsilon)^2/sigmastar1 - sigx9_1/sigma_sq1))/sigma_sq1 +
      sigx2_1 * sigmastar1/ewu1 - sigx17_1^2 * pwZ * sigx1_1/sigx19) *
    pwZ * sigx1_1/sigx19, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx17_1 * sigx17_2 * (1 - pwZ) * pwZ *
      sigx1_1 * sigx1_2/sigx19^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (vsq1 * dmusig1) +
      0.5 * ((dmusig1 * ewv/sigma_sq1 - S^2 * vsq1 * sigx6_1 *
        dmusig1 * ewu1 * (epsilon)^2/sigmastar1) * ewu1/sigma_sq1 -
        0.5 * (usq1 * vsq1 * dmusig1)))/sigmastar1 +
      dpsv1 * usq1 * sigx14_1 + (S * (pmusig1 - ewu1 *
      (pmusig1/sigma_sq1 + S * sigx6_1 * dmusig1 * (epsilon))) *
      (epsilon) - sigx11_1 * ewu1)/sigma_sq1)/sigma_sq1 -
      S^2 * (((0.5 * (usq1 * ewv/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
        1) * ewv/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) *
        ewu1/sigmastar1 + 2 * sigx12_1)/sigx13_1 - ((starsq1^2 +
        2 * (sigx5_1 * sigma_sq1^2 * sigmastar1)) * sigmastar1 +
        0.5 * (starsq1^2 * usq1 * ewv/sigmastar1)) *
        sigx12_1 * ewu1/sigx13_1^2) * sigx2_1 * ewu1 *
        (epsilon)^2) * sigmastar1 + (0.5 * (sigx15_1 *
      usq1 * ewv) + S^2 * sigx16_1 * sigx6_1 * ewu1 * (epsilon)^2)/starsq1 +
      0.5 * (((dpsv1 * usq1 * vsq1 + sigx2_1 * ewv/sigma_sq1) *
        ewu1/sigma_sq1 + vsq1 * sigx2_1)/starsq1 - sigx5_1 *
        vsq1 * sigx2_1 * ewu1/starsq1^2) - (sigx20 *
      sigx17_1/sigx19 + sigx16_1/ewu1)) * pwZ * sigx1_1 *
      ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx17_1 * (1 - sigx21 *
      pwZ/sigx19) * dwZ * sigx1_1/sigx19, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((dpsv2 * usq2 +
      S * ewu2 * pmusig2 * (epsilon)/sigma_sq2 - dmusig2 *
      sigmastar2) * ewu2/sigma_sq2 - 0.5 * (usq2 * sigx2_2)) +
      0.5 * (sigx7_2 * ewu2/sigma_sq2) - 0.5 * sigx2_2) *
      ewv/sigmastar2 - dpsv2 * sigmastar2) * usq2 + (usq2 *
      (ewu2 * (S^2 * sigx6_2 * dmusig2 * (epsilon)^2 -
        dpsv2/sigma_sq2) - 0.5 * (0.5 * (usq2 * dmusig2 *
        ewv/sigmastar2) + S^2 * sigx6_2 * dmusig2 * ewu2 *
        (epsilon)^2)) + S^2 * ((dpsv2 * usq2 * sigx6_2/sigma_sq2 -
      ((0.5 * (ewu2/sigma_sq2) + 1 - 0.5 * (0.5 * usq2 +
        ewu2/sigma_sq2)) * usq2 * ewv/sigmastar2 + (2 -
        2 * (sigx5_2^2 * ewu2 * sigma_sq2/starsq2^2)) *
        sigmastar2) * sigx2_2/starsq2^2) * ewu2 + (1 -
      0.5 * usq2) * sigx6_2 * sigx2_2) * ewu2 * (epsilon)^2/sigmastar2) *
      sigmastar2 + ewu2 * (S^2 * sigx17_2 * sigx6_2 * ewu2 *
      (epsilon)^2/sigmastar2 - sigx9_2/sigma_sq2))/sigma_sq2 +
      sigx2_2 * sigmastar2/ewu2 - sigx17_2^2 * (1 - pwZ) *
      sigx1_2/sigx19) * (1 - pwZ) * sigx1_2/sigx19, FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (vsq2 * dmusig2) +
      0.5 * ((dmusig2 * ewv/sigma_sq2 - S^2 * vsq2 * sigx6_2 *
        dmusig2 * ewu2 * (epsilon)^2/sigmastar2) * ewu2/sigma_sq2 -
        0.5 * (usq2 * vsq2 * dmusig2)))/sigmastar2 +
      dpsv2 * usq2 * sigx14_2 + (S * (pmusig2 - ewu2 *
      (pmusig2/sigma_sq2 + S * sigx6_2 * dmusig2 * (epsilon))) *
      (epsilon) - sigx11_2 * ewu2)/sigma_sq2)/sigma_sq2 -
      S^2 * (((0.5 * (usq2 * ewv/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
        1) * ewv/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
        ewu2/sigmastar2 + 2 * sigx12_2)/sigx13_2 - ((starsq2^2 +
        2 * (sigx5_2 * sigma_sq2^2 * sigmastar2)) * sigmastar2 +
        0.5 * (starsq2^2 * usq2 * ewv/sigmastar2)) *
        sigx12_2 * ewu2/sigx13_2^2) * sigx2_2 * ewu2 *
        (epsilon)^2) * sigmastar2 + (0.5 * (sigx15_2 *
      usq2 * ewv) + S^2 * sigx16_2 * sigx6_2 * ewu2 * (epsilon)^2)/starsq2 +
      0.5 * (((dpsv2 * usq2 * vsq2 + sigx2_2 * ewv/sigma_sq2) *
        ewu2/sigma_sq2 + vsq2 * sigx2_2)/starsq2 - sigx5_2 *
        vsq2 * sigx2_2 * ewu2/starsq2^2) - (sigx20 *
      sigx17_2/sigx19 + sigx16_2/ewu2)) * (1 - pwZ) * sigx1_2 *
      ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx17_2 * (sigx21 * (1 - pwZ)/sigx19 + 1) * dwZ * sigx1_2/sigx19),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((sigx11_1 * sigx14_1 +
      (S * (S * sigx12_1 * dmusig1 * ewu1 * (epsilon)/starsq1^2 -
        2 * (pmusig1/sigma_sq1)) * (epsilon) - 0.5 *
        sigx10_1)/sigma_sq1) * ewv + (0.5 * (ewv * (S^2 *
      sigx12_1 * dmusig1 * ewu1^2 * (epsilon)^2/(starsq1^2 *
      sigmastar1) - dmusig1)/sigma_sq1 - 0.5 * (vsq1 *
      dmusig1)) + 0.5 * dmusig1) * vsq1/sigmastar1 + S *
      pmusig1 * (epsilon)/sigma_sq1) * ewu1/sigma_sq1 +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx12_1 * (1/sigx13_1 - ((starsq1^2 +
        2 * (sigx12_1 * sigma_sq1^2 * sigmastar1)) *
        sigmastar1 + 0.5 * (starsq1^2 * vsq1 * ewu1/sigmastar1)) *
        ewv/sigx13_1^2) + (0.5 * (ewv/sigma_sq1) - 0.5 *
        (0.5 * vsq1 + ewv/sigma_sq1)) * vsq1/(starsq1^2 *
        ewv)) * ewu1^2 * (epsilon)^2) * sigx2_1) * sigmastar1 +
      (sigx16_1 * sigx14_1 + (0.5 * (((0.5 * sigx10_1 +
        2 * (S * pmusig1 * (epsilon)/sigma_sq1)) * ewu1 -
        dmusig1 * sigmastar1)/(sigma_sq1^2 * sigmastar1) -
        sigx12_1 * sigx2_1/starsq1^2) + 0.5 * (sigx15_1/starsq1)) *
        vsq1 * ewu1) * ewv + 0.5 * (vsq1 * sigx2_1 *
      ewu1/starsq1)) * pwZ * sigx1_1/ewu1 + ((((sigx11_2 *
      sigx14_2 + (S * (S * sigx12_2 * dmusig2 * ewu2 *
      (epsilon)/starsq2^2 - 2 * (pmusig2/sigma_sq2)) *
      (epsilon) - 0.5 * sigx10_2)/sigma_sq2) * ewv + (0.5 *
      (ewv * (S^2 * sigx12_2 * dmusig2 * ewu2^2 * (epsilon)^2/(starsq2^2 *
        sigmastar2) - dmusig2)/sigma_sq2 - 0.5 * (vsq2 *
        dmusig2)) + 0.5 * dmusig2) * vsq2/sigmastar2 +
      S * pmusig2 * (epsilon)/sigma_sq2) * ewu2/sigma_sq2 +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx12_2 * (1/sigx13_2 - ((starsq2^2 +
        2 * (sigx12_2 * sigma_sq2^2 * sigmastar2)) *
        sigmastar2 + 0.5 * (starsq2^2 * vsq2 * ewu2/sigmastar2)) *
        ewv/sigx13_2^2) + (0.5 * (ewv/sigma_sq2) - 0.5 *
        (0.5 * vsq2 + ewv/sigma_sq2)) * vsq2/(starsq2^2 *
        ewv)) * ewu2^2 * (epsilon)^2) * sigx2_2) * sigmastar2 +
      (sigx16_2 * sigx14_2 + (0.5 * (((0.5 * sigx10_2 +
        2 * (S * pmusig2 * (epsilon)/sigma_sq2)) * ewu2 -
        dmusig2 * sigmastar2)/(sigma_sq2^2 * sigmastar2) -
        sigx12_2 * sigx2_2/starsq2^2) + 0.5 * (sigx15_2/starsq2)) *
        vsq2 * ewu2) * ewv + 0.5 * (vsq2 * sigx2_2 *
      ewu2/starsq2)) * (1 - pwZ) * sigx1_2/ewu2 - sigx20^2 *
      ewv/sigx19) * ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (sigx16_1 * sigx1_1/ewu1 - (sigx20 *
      sigx21/sigx19 + sigx16_2 * sigx1_2/ewu2)) * dwZ *
      ewv/sigx19, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx21 * dwZ/sigx19 +
      Wz) * sigx21 * dwZ/sigx19), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cloglog specification class membership
chessmisfraynormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv)
  sigma_sq2 <- (ewu2 + ewv)
  sigmastar1 <- sqrt(ewu1 * ewv/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  vsq1 <- (1 - ewv/sigma_sq1)
  vsq2 <- (1 - ewv/sigma_sq2)
  musig1 <- (S * ewu1 * (epsilon)/starsq1)
  musig2 <- (S * ewu2 * (epsilon)/starsq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  wuv1 <- (sigmastar1/ewv - ewu1/starsq1)
  wuv2 <- (sigmastar2/ewv - ewu2/starsq2)
  dvs1 <- (dmusig1 * ewv/sigmastar1)
  dvs2 <- (dmusig2 * ewv/sigmastar2)
  dpsv1 <- (0.5 * dvs1 - S * pmusig1 * (epsilon))
  dpsv2 <- (0.5 * dvs2 - S * pmusig2 * (epsilon))
  vepsi <- ((S * (epsilon))^2/(2 * ewv)^2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv))
  sigx2_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/sigma_sq1)
  sigx2_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/sigma_sq2)
  sigx3_1 <- (pmusig1 + S * dmusig1 * wuv1 * (epsilon))
  sigx3_2 <- (pmusig2 + S * dmusig2 * wuv2 * (epsilon))
  sigx4_1 <- (ewu1 * sigx3_1/sigma_sq1 + S * usq1 * sigx2_1 *
    (epsilon)/ewv)
  sigx4_2 <- (ewu2 * sigx3_2/sigma_sq2 + S * usq2 * sigx2_2 *
    (epsilon)/ewv)
  sigx5_1 <- (0.5 * (usq1 * ewv/sigmastar1) + sigmastar1)
  sigx5_2 <- (0.5 * (usq2 * ewv/sigmastar2) + sigmastar2)
  sigx6_1 <- (1/starsq1 - sigx5_1 * ewu1/starsq1^2)
  sigx6_2 <- (1/starsq2 - sigx5_2 * ewu2/starsq2^2)
  sigx7_1 <- (dpsv1 * usq1 + S^2 * sigx6_1 * sigx2_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx7_2 <- (dpsv2 * usq2 + S^2 * sigx6_2 * sigx2_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx8_1 <- (usq1 * sigx2_1 * ewv/sigmastar1)
  sigx8_2 <- (usq2 * sigx2_2 * ewv/sigmastar2)
  sigx9_1 <- (sigx7_1 * sigmastar1 + 0.5 * sigx8_1)
  sigx9_2 <- (sigx7_2 * sigmastar2 + 0.5 * sigx8_2)
  sigx10_1 <- (vsq1 * dmusig1/sigmastar1)
  sigx10_2 <- (vsq2 * dmusig2/sigmastar2)
  sigx11_1 <- (0.5 * sigx10_1 + S * pmusig1 * (epsilon)/sigma_sq1)
  sigx11_2 <- (0.5 * sigx10_2 + S * pmusig2 * (epsilon)/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (starsq1^2 * sigma_sq1 * sigmastar1)
  sigx13_2 <- (starsq2^2 * sigma_sq2 * sigmastar2)
  sigx14_1 <- (2 * vepsi - S^2 * sigx12_1 * ewu1^2 * (epsilon)^2/sigx13_1)
  sigx14_2 <- (2 * vepsi - S^2 * sigx12_2 * ewu2^2 * (epsilon)^2/sigx13_2)
  sigx15_1 <- (sigx11_1 * ewu1/sigma_sq1 + sigx14_1 * sigx2_1)
  sigx15_2 <- (sigx11_2 * ewu2/sigma_sq2 + sigx14_2 * sigx2_2)
  sigx16_1 <- (sigx15_1 * sigmastar1 + 0.5 * (vsq1 * sigx2_1 *
    ewu1/starsq1))
  sigx16_2 <- (sigx15_2 * sigmastar2 + 0.5 * (vsq2 * sigx2_2 *
    ewu2/starsq2))
  sigx17_1 <- (sigx9_1/sigma_sq1 - sigx2_1 * sigmastar1/ewu1)
  sigx17_2 <- (sigx9_2/sigma_sq2 - sigx2_2 * sigmastar2/ewu2)
  sigx18 <- ((1 - prZ) * sigx1_1 * sigx4_1 * sigmastar1/ewu1 +
    prZ * sigx1_2 * sigx4_2 * sigmastar2/ewu2)
  sigx19 <- ((1 - prZ) * sigx2_1 * sigx1_1 * sigmastar1/ewu1 +
    sigx2_2 * prZ * sigx1_2 * sigmastar2/ewu2)
  sigx20 <- (sigx16_1 * (1 - prZ) * sigx1_1/ewu1 + sigx16_2 *
    prZ * sigx1_2/ewu2)
  sigx21 <- (sigx2_1 * sigx1_1 * sigmastar1/ewu1 - sigx2_2 *
    sigx1_2 * sigmastar2/ewu2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((((S^2 * ewu1 * (epsilon)^2/(sigma_sq1 *
      ewv) - 1) * wuv1 + ewu1/starsq1) * dmusig1 * ewu1/sigma_sq1 +
      usq1 * (S * ((2 * (S * dmusig1 * wuv1 * (epsilon)) +
        3 * pmusig1) * ewu1/sigma_sq1 + S * usq1 * sigx2_1 *
        (epsilon)/ewv) * (epsilon) - dmusig1 * sigmastar1)/ewv) *
      (1 - prZ) * sigx1_1 * sigmastar1/ewu1 + (((S^2 *
      ewu2 * (epsilon)^2/(sigma_sq2 * ewv) - 1) * wuv2 +
      ewu2/starsq2) * dmusig2 * ewu2/sigma_sq2 + usq2 *
      (S * ((2 * (S * dmusig2 * wuv2 * (epsilon)) + 3 *
        pmusig2) * ewu2/sigma_sq2 + S * usq2 * sigx2_2 *
        (epsilon)/ewv) * (epsilon) - dmusig2 * sigmastar2)/ewv) *
      prZ * sigx1_2 * sigmastar2/ewu2 - sigx18^2/sigx19)/sigx19,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx17_1 * (S * usq1 *
      (epsilon)/ewv - sigx18/sigx19) + ((usq1 * (pmusig1 -
      0.5 * (S * dmusig1 * ewu1 * (epsilon)/starsq1)) +
      S * sigx6_1 * ewu1 * (S * ewu1 * sigx3_1 * (epsilon)/sigma_sq1 -
        2 * sigx2_1) * (epsilon)/sigmastar1) * sigmastar1 +
      (0.5 * (usq1 * ewu1 * ewv/starsq1) - sigmastar1) *
        sigx3_1)/sigma_sq1) * (1 - prZ) * sigx1_1/sigx19,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx17_2 * (S * usq2 *
      (epsilon)/ewv - sigx18/sigx19) + ((usq2 * (pmusig2 -
      0.5 * (S * dmusig2 * ewu2 * (epsilon)/starsq2)) +
      S * sigx6_2 * ewu2 * (S * ewu2 * sigx3_2 * (epsilon)/sigma_sq2 -
        2 * sigx2_2) * (epsilon)/sigmastar2) * sigmastar2 +
      (0.5 * (usq2 * ewu2 * ewv/starsq2) - sigmastar2) *
        sigx3_2)/sigma_sq2) * prZ * sigx1_2/sigx19, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((sigx14_1 * sigx3_1 + (S * (0.5 * (vsq1/ewv) +
    1/sigma_sq1) * dmusig1 * ewu1 * (epsilon)/sigmastar1 -
    pmusig1)/sigma_sq1) * ewu1/sigma_sq1 + S * (2 * (sigx12_1 *
    ewu1^2/sigx13_1) - 4/(2 * ewv)^2) * sigx2_1 * (epsilon)) *
    sigmastar1 + 0.5 * (vsq1 * ewu1^2 * sigx3_1/(sigma_sq1^2 *
    sigmastar1)) + S * sigx16_1 * usq1 * (epsilon)/ewv) *
    (1 - prZ) * sigx1_1/ewu1 + (((sigx14_2 * sigx3_2 + (S *
    (0.5 * (vsq2/ewv) + 1/sigma_sq2) * dmusig2 * ewu2 * (epsilon)/sigmastar2 -
    pmusig2)/sigma_sq2) * ewu2/sigma_sq2 + S * (2 * (sigx12_2 *
    ewu2^2/sigx13_2) - 4/(2 * ewv)^2) * sigx2_2 * (epsilon)) *
    sigmastar2 + 0.5 * (vsq2 * ewu2^2 * sigx3_2/(sigma_sq2^2 *
    sigmastar2)) + S * sigx16_2 * usq2 * (epsilon)/ewv) *
    prZ * sigx1_2/ewu2 - sigx20 * sigx18/sigx19) * ewv/sigx19,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * prZ * (sigx1_1 * sigx4_1 *
      sigmastar1/ewu1 - (sigx18 * sigx21/sigx19 + sigx1_2 *
      sigx4_2 * sigmastar2/ewu2)) * ewz/sigx19, FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((dpsv1 * usq1 + S * ewu1 * pmusig1 * (epsilon)/sigma_sq1 -
      dmusig1 * sigmastar1) * ewu1/sigma_sq1 - 0.5 * (usq1 *
      sigx2_1)) + 0.5 * (sigx7_1 * ewu1/sigma_sq1) - 0.5 *
      sigx2_1) * ewv/sigmastar1 - dpsv1 * sigmastar1) *
      usq1 + (usq1 * (ewu1 * (S^2 * sigx6_1 * dmusig1 *
      (epsilon)^2 - dpsv1/sigma_sq1) - 0.5 * (0.5 * (usq1 *
      dmusig1 * ewv/sigmastar1) + S^2 * sigx6_1 * dmusig1 *
      ewu1 * (epsilon)^2)) + S^2 * ((dpsv1 * usq1 * sigx6_1/sigma_sq1 -
      ((0.5 * (ewu1/sigma_sq1) + 1 - 0.5 * (0.5 * usq1 +
        ewu1/sigma_sq1)) * usq1 * ewv/sigmastar1 + (2 -
        2 * (sigx5_1^2 * ewu1 * sigma_sq1/starsq1^2)) *
        sigmastar1) * sigx2_1/starsq1^2) * ewu1 + (1 -
      0.5 * usq1) * sigx6_1 * sigx2_1) * ewu1 * (epsilon)^2/sigmastar1) *
      sigmastar1 + ewu1 * (S^2 * sigx17_1 * sigx6_1 * ewu1 *
      (epsilon)^2/sigmastar1 - sigx9_1/sigma_sq1))/sigma_sq1 +
      sigx2_1 * sigmastar1/ewu1 - sigx17_1^2 * (1 - prZ) *
      sigx1_1/sigx19) * (1 - prZ) * sigx1_1/sigx19, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx17_1 * sigx17_2 * prZ * (1 - prZ) *
      sigx1_1 * sigx1_2/sigx19^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (vsq1 * dmusig1) +
      0.5 * ((dmusig1 * ewv/sigma_sq1 - S^2 * vsq1 * sigx6_1 *
        dmusig1 * ewu1 * (epsilon)^2/sigmastar1) * ewu1/sigma_sq1 -
        0.5 * (usq1 * vsq1 * dmusig1)))/sigmastar1 +
      dpsv1 * usq1 * sigx14_1 + (S * (pmusig1 - ewu1 *
      (pmusig1/sigma_sq1 + S * sigx6_1 * dmusig1 * (epsilon))) *
      (epsilon) - sigx11_1 * ewu1)/sigma_sq1)/sigma_sq1 -
      S^2 * (((0.5 * (usq1 * ewv/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 -
        1) * ewv/sigma_sq1 + 1 - 0.5 * (usq1 * vsq1))) *
        ewu1/sigmastar1 + 2 * sigx12_1)/sigx13_1 - ((starsq1^2 +
        2 * (sigx5_1 * sigma_sq1^2 * sigmastar1)) * sigmastar1 +
        0.5 * (starsq1^2 * usq1 * ewv/sigmastar1)) *
        sigx12_1 * ewu1/sigx13_1^2) * sigx2_1 * ewu1 *
        (epsilon)^2) * sigmastar1 + (0.5 * (sigx15_1 *
      usq1 * ewv) + S^2 * sigx16_1 * sigx6_1 * ewu1 * (epsilon)^2)/starsq1 +
      0.5 * (((dpsv1 * usq1 * vsq1 + sigx2_1 * ewv/sigma_sq1) *
        ewu1/sigma_sq1 + vsq1 * sigx2_1)/starsq1 - sigx5_1 *
        vsq1 * sigx2_1 * ewu1/starsq1^2) - (sigx20 *
      sigx17_1/sigx19 + sigx16_1/ewu1)) * (1 - prZ) * sigx1_1 *
      ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx17_1 * (1 - sigx21 *
      (1 - prZ)/sigx19) * prZ * sigx1_1 * ewz/sigx19, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((dpsv2 * usq2 +
      S * ewu2 * pmusig2 * (epsilon)/sigma_sq2 - dmusig2 *
      sigmastar2) * ewu2/sigma_sq2 - 0.5 * (usq2 * sigx2_2)) +
      0.5 * (sigx7_2 * ewu2/sigma_sq2) - 0.5 * sigx2_2) *
      ewv/sigmastar2 - dpsv2 * sigmastar2) * usq2 + (usq2 *
      (ewu2 * (S^2 * sigx6_2 * dmusig2 * (epsilon)^2 -
        dpsv2/sigma_sq2) - 0.5 * (0.5 * (usq2 * dmusig2 *
        ewv/sigmastar2) + S^2 * sigx6_2 * dmusig2 * ewu2 *
        (epsilon)^2)) + S^2 * ((dpsv2 * usq2 * sigx6_2/sigma_sq2 -
      ((0.5 * (ewu2/sigma_sq2) + 1 - 0.5 * (0.5 * usq2 +
        ewu2/sigma_sq2)) * usq2 * ewv/sigmastar2 + (2 -
        2 * (sigx5_2^2 * ewu2 * sigma_sq2/starsq2^2)) *
        sigmastar2) * sigx2_2/starsq2^2) * ewu2 + (1 -
      0.5 * usq2) * sigx6_2 * sigx2_2) * ewu2 * (epsilon)^2/sigmastar2) *
      sigmastar2 + ewu2 * (S^2 * sigx17_2 * sigx6_2 * ewu2 *
      (epsilon)^2/sigmastar2 - sigx9_2/sigma_sq2))/sigma_sq2 +
      sigx2_2 * sigmastar2/ewu2 - sigx17_2^2 * prZ * sigx1_2/sigx19) *
      prZ * sigx1_2/sigx19, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (vsq2 * dmusig2) +
      0.5 * ((dmusig2 * ewv/sigma_sq2 - S^2 * vsq2 * sigx6_2 *
        dmusig2 * ewu2 * (epsilon)^2/sigmastar2) * ewu2/sigma_sq2 -
        0.5 * (usq2 * vsq2 * dmusig2)))/sigmastar2 +
      dpsv2 * usq2 * sigx14_2 + (S * (pmusig2 - ewu2 *
      (pmusig2/sigma_sq2 + S * sigx6_2 * dmusig2 * (epsilon))) *
      (epsilon) - sigx11_2 * ewu2)/sigma_sq2)/sigma_sq2 -
      S^2 * (((0.5 * (usq2 * ewv/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
        1) * ewv/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) *
        ewu2/sigmastar2 + 2 * sigx12_2)/sigx13_2 - ((starsq2^2 +
        2 * (sigx5_2 * sigma_sq2^2 * sigmastar2)) * sigmastar2 +
        0.5 * (starsq2^2 * usq2 * ewv/sigmastar2)) *
        sigx12_2 * ewu2/sigx13_2^2) * sigx2_2 * ewu2 *
        (epsilon)^2) * sigmastar2 + (0.5 * (sigx15_2 *
      usq2 * ewv) + S^2 * sigx16_2 * sigx6_2 * ewu2 * (epsilon)^2)/starsq2 +
      0.5 * (((dpsv2 * usq2 * vsq2 + sigx2_2 * ewv/sigma_sq2) *
        ewu2/sigma_sq2 + vsq2 * sigx2_2)/starsq2 - sigx5_2 *
        vsq2 * sigx2_2 * ewu2/starsq2^2) - (sigx20 *
      sigx17_2/sigx19 + sigx16_2/ewu2)) * prZ * sigx1_2 *
      ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx17_2 * (sigx21 * prZ/sigx19 + 1) * prZ * sigx1_2 *
      ewz/sigx19), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((sigx11_1 * sigx14_1 +
      (S * (S * sigx12_1 * dmusig1 * ewu1 * (epsilon)/starsq1^2 -
        2 * (pmusig1/sigma_sq1)) * (epsilon) - 0.5 *
        sigx10_1)/sigma_sq1) * ewv + (0.5 * (ewv * (S^2 *
      sigx12_1 * dmusig1 * ewu1^2 * (epsilon)^2/(starsq1^2 *
      sigmastar1) - dmusig1)/sigma_sq1 - 0.5 * (vsq1 *
      dmusig1)) + 0.5 * dmusig1) * vsq1/sigmastar1 + S *
      pmusig1 * (epsilon)/sigma_sq1) * ewu1/sigma_sq1 +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx12_1 * (1/sigx13_1 - ((starsq1^2 +
        2 * (sigx12_1 * sigma_sq1^2 * sigmastar1)) *
        sigmastar1 + 0.5 * (starsq1^2 * vsq1 * ewu1/sigmastar1)) *
        ewv/sigx13_1^2) + (0.5 * (ewv/sigma_sq1) - 0.5 *
        (0.5 * vsq1 + ewv/sigma_sq1)) * vsq1/(starsq1^2 *
        ewv)) * ewu1^2 * (epsilon)^2) * sigx2_1) * sigmastar1 +
      (sigx16_1 * sigx14_1 + (0.5 * (((0.5 * sigx10_1 +
        2 * (S * pmusig1 * (epsilon)/sigma_sq1)) * ewu1 -
        dmusig1 * sigmastar1)/(sigma_sq1^2 * sigmastar1) -
        sigx12_1 * sigx2_1/starsq1^2) + 0.5 * (sigx15_1/starsq1)) *
        vsq1 * ewu1) * ewv + 0.5 * (vsq1 * sigx2_1 *
      ewu1/starsq1)) * (1 - prZ) * sigx1_1/ewu1 + ((((sigx11_2 *
      sigx14_2 + (S * (S * sigx12_2 * dmusig2 * ewu2 *
      (epsilon)/starsq2^2 - 2 * (pmusig2/sigma_sq2)) *
      (epsilon) - 0.5 * sigx10_2)/sigma_sq2) * ewv + (0.5 *
      (ewv * (S^2 * sigx12_2 * dmusig2 * ewu2^2 * (epsilon)^2/(starsq2^2 *
        sigmastar2) - dmusig2)/sigma_sq2 - 0.5 * (vsq2 *
        dmusig2)) + 0.5 * dmusig2) * vsq2/sigmastar2 +
      S * pmusig2 * (epsilon)/sigma_sq2) * ewu2/sigma_sq2 +
      ((2 - 16 * (ewv^2/(2 * ewv)^2)) * (S * (epsilon))^2/(2 *
        ewv)^2 - S^2 * (sigx12_2 * (1/sigx13_2 - ((starsq2^2 +
        2 * (sigx12_2 * sigma_sq2^2 * sigmastar2)) *
        sigmastar2 + 0.5 * (starsq2^2 * vsq2 * ewu2/sigmastar2)) *
        ewv/sigx13_2^2) + (0.5 * (ewv/sigma_sq2) - 0.5 *
        (0.5 * vsq2 + ewv/sigma_sq2)) * vsq2/(starsq2^2 *
        ewv)) * ewu2^2 * (epsilon)^2) * sigx2_2) * sigmastar2 +
      (sigx16_2 * sigx14_2 + (0.5 * (((0.5 * sigx10_2 +
        2 * (S * pmusig2 * (epsilon)/sigma_sq2)) * ewu2 -
        dmusig2 * sigmastar2)/(sigma_sq2^2 * sigmastar2) -
        sigx12_2 * sigx2_2/starsq2^2) + 0.5 * (sigx15_2/starsq2)) *
        vsq2 * ewu2) * ewv + 0.5 * (vsq2 * sigx2_2 *
      ewu2/starsq2)) * prZ * sigx1_2/ewu2 - sigx20^2 *
      ewv/sigx19) * ewv/sigx19, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (sigx16_1 * sigx1_1/ewu1 - (sigx20 *
      sigx21/sigx19 + sigx16_2 * sigx1_2/ewu2)) * prZ *
      ewv * ewz/sigx19, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * sigx21 * (1 - (sigx21 * prZ/sigx19 +
      1) * ewz) * prZ * ewz/sigx19, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}


# Optimization using different algorithms ----------
#' optimizations solve for misf rayleigh-normal distribution
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
misfraynormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfraynormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmisfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfraynormlike_logit,
      grad = cgradmisfraynormlike_logit, hess = chessmisfraynormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfraynormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfraynormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfraynormlike_logit(mleObj$par,
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
      mleObj$hessian <- chessmisfraynormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfraynormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfraynormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfraynormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initRay = initRay))
}

# cauchit specification class membership
misfraynormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfraynormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmisfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfraynormlike_cauchit,
      grad = cgradmisfraynormlike_cauchit, hess = chessmisfraynormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfraynormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfraynormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfraynormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- chessmisfraynormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfraynormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfraynormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfraynormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initRay = initRay))
}

# probit specification class membership
misfraynormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfraynormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmisfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfraynormlike_probit,
      grad = cgradmisfraynormlike_probit, hess = chessmisfraynormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfraynormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfraynormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfraynormlike_probit(mleObj$par,
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
      mleObj$hessian <- chessmisfraynormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfraynormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfraynormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfraynormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initRay = initRay))
}

# cloglog specification class membership
misfraynormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfraynormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmisfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfraynormlike_cloglog,
      grad = cgradmisfraynormlike_cloglog, hess = chessmisfraynormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfraynormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfraynormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfraynormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- chessmisfraynormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfraynormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfraynormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfraynormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initRay = initRay))
}

# Conditional efficiencies estimation ----------
#' efficiencies for misf rayleigh-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# logit specification class membership
cmisfraynormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- (mustar1 * sigmastar1 * dnorm(mustar1/sigmastar1) +
    (mustar1^2 + sigmastar1^2) * pnorm(mustar1/sigmastar1))/(sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  u_c2 <- (mustar2 * sigmastar2 * dnorm(mustar2/sigmastar2) +
    (mustar2^2 + sigmastar2^2) * pnorm(mustar2/sigmastar2))/(sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + sigmastar1^2/2) * (sigmastar1 *
      dnorm(mustar1/sigmastar1 - sigmastar1) + (mustar1 -
      sigmastar1^2) * pnorm(mustar1/sigmastar1 - sigmastar1))/(sigmastar1 *
      dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
    teBC_c2 <- exp(-mustar2 + sigmastar2^2/2) * (sigmastar2 *
      dnorm(mustar2/sigmastar2 - sigmastar2) + (mustar2 -
      sigmastar2^2) * pnorm(mustar2/sigmastar2 - sigmastar2))/(sigmastar2 *
      dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + sigmastar1^2/2) *
      (sigmastar1 * dnorm(mustar1/sigmastar1 + sigmastar1) +
        (mustar1 + sigmastar1^2) * pnorm(mustar1/sigmastar1 +
          sigmastar1))/(sigmastar1 * dnorm(mustar1/sigmastar1) +
      mustar1 * pnorm(mustar1/sigmastar1))
    teBC_reciprocal_c2 <- exp(mustar2 + sigmastar2^2/2) *
      (sigmastar2 * dnorm(mustar2/sigmastar2 + sigmastar2) +
        (mustar2 + sigmastar2^2) * pnorm(mustar2/sigmastar2 +
          sigmastar2))/(sigmastar2 * dnorm(mustar2/sigmastar2) +
      mustar2 * pnorm(mustar2/sigmastar2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# cauchit specification class membership
cmisfraynormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- (mustar1 * sigmastar1 * dnorm(mustar1/sigmastar1) +
    (mustar1^2 + sigmastar1^2) * pnorm(mustar1/sigmastar1))/(sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  u_c2 <- (mustar2 * sigmastar2 * dnorm(mustar2/sigmastar2) +
    (mustar2^2 + sigmastar2^2) * pnorm(mustar2/sigmastar2))/(sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + sigmastar1^2/2) * (sigmastar1 *
      dnorm(mustar1/sigmastar1 - sigmastar1) + (mustar1 -
      sigmastar1^2) * pnorm(mustar1/sigmastar1 - sigmastar1))/(sigmastar1 *
      dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
    teBC_c2 <- exp(-mustar2 + sigmastar2^2/2) * (sigmastar2 *
      dnorm(mustar2/sigmastar2 - sigmastar2) + (mustar2 -
      sigmastar2^2) * pnorm(mustar2/sigmastar2 - sigmastar2))/(sigmastar2 *
      dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + sigmastar1^2/2) *
      (sigmastar1 * dnorm(mustar1/sigmastar1 + sigmastar1) +
        (mustar1 + sigmastar1^2) * pnorm(mustar1/sigmastar1 +
          sigmastar1))/(sigmastar1 * dnorm(mustar1/sigmastar1) +
      mustar1 * pnorm(mustar1/sigmastar1))
    teBC_reciprocal_c2 <- exp(mustar2 + sigmastar2^2/2) *
      (sigmastar2 * dnorm(mustar2/sigmastar2 + sigmastar2) +
        (mustar2 + sigmastar2^2) * pnorm(mustar2/sigmastar2 +
          sigmastar2))/(sigmastar2 * dnorm(mustar2/sigmastar2) +
      mustar2 * pnorm(mustar2/sigmastar2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# probit specification class membership
cmisfraynormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- (mustar1 * sigmastar1 * dnorm(mustar1/sigmastar1) +
    (mustar1^2 + sigmastar1^2) * pnorm(mustar1/sigmastar1))/(sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  u_c2 <- (mustar2 * sigmastar2 * dnorm(mustar2/sigmastar2) +
    (mustar2^2 + sigmastar2^2) * pnorm(mustar2/sigmastar2))/(sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + sigmastar1^2/2) * (sigmastar1 *
      dnorm(mustar1/sigmastar1 - sigmastar1) + (mustar1 -
      sigmastar1^2) * pnorm(mustar1/sigmastar1 - sigmastar1))/(sigmastar1 *
      dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
    teBC_c2 <- exp(-mustar2 + sigmastar2^2/2) * (sigmastar2 *
      dnorm(mustar2/sigmastar2 - sigmastar2) + (mustar2 -
      sigmastar2^2) * pnorm(mustar2/sigmastar2 - sigmastar2))/(sigmastar2 *
      dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + sigmastar1^2/2) *
      (sigmastar1 * dnorm(mustar1/sigmastar1 + sigmastar1) +
        (mustar1 + sigmastar1^2) * pnorm(mustar1/sigmastar1 +
          sigmastar1))/(sigmastar1 * dnorm(mustar1/sigmastar1) +
      mustar1 * pnorm(mustar1/sigmastar1))
    teBC_reciprocal_c2 <- exp(mustar2 + sigmastar2^2/2) *
      (sigmastar2 * dnorm(mustar2/sigmastar2 + sigmastar2) +
        (mustar2 + sigmastar2^2) * pnorm(mustar2/sigmastar2 +
          sigmastar2))/(sigmastar2 * dnorm(mustar2/sigmastar2) +
      mustar2 * pnorm(mustar2/sigmastar2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# cloglog specification class membership
cmisfraynormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- (mustar1 * sigmastar1 * dnorm(mustar1/sigmastar1) +
    (mustar1^2 + sigmastar1^2) * pnorm(mustar1/sigmastar1))/(sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  u_c2 <- (mustar2 * sigmastar2 * dnorm(mustar2/sigmastar2) +
    (mustar2^2 + sigmastar2^2) * pnorm(mustar2/sigmastar2))/(sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + sigmastar1^2/2) * (sigmastar1 *
      dnorm(mustar1/sigmastar1 - sigmastar1) + (mustar1 -
      sigmastar1^2) * pnorm(mustar1/sigmastar1 - sigmastar1))/(sigmastar1 *
      dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
    teBC_c2 <- exp(-mustar2 + sigmastar2^2/2) * (sigmastar2 *
      dnorm(mustar2/sigmastar2 - sigmastar2) + (mustar2 -
      sigmastar2^2) * pnorm(mustar2/sigmastar2 - sigmastar2))/(sigmastar2 *
      dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + sigmastar1^2/2) *
      (sigmastar1 * dnorm(mustar1/sigmastar1 + sigmastar1) +
        (mustar1 + sigmastar1^2) * pnorm(mustar1/sigmastar1 +
          sigmastar1))/(sigmastar1 * dnorm(mustar1/sigmastar1) +
      mustar1 * pnorm(mustar1/sigmastar1))
    teBC_reciprocal_c2 <- exp(mustar2 + sigmastar2^2/2) *
      (sigmastar2 * dnorm(mustar2/sigmastar2 + sigmastar2) +
        (mustar2 + sigmastar2^2) * pnorm(mustar2/sigmastar2 +
          sigmastar2))/(sigmastar2 * dnorm(mustar2/sigmastar2) +
      mustar2 * pnorm(mustar2/sigmastar2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for misf rayleigh-normal distribution
#' @param object object of class sfacross
#' @noRd
# logit specification class membership
cmisfmargraynorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargraynorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (4 - pi)/2, ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (4 - pi)/2, ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cauchit specification class membership
cmisfmargraynorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargraynorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (4 - pi)/2, ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (4 - pi)/2, ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# probit specification class membership
cmisfmargraynorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargraynorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (4 - pi)/2, ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (4 - pi)/2, ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cloglog specification class membership
cmisfmargraynorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmargraynorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
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
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv)/(exp(Wu1) + exp(Wv)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv)/(exp(Wu2) + exp(Wv)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv)))/(exp(Wv/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (4 - pi)/2, ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (4 - pi)/2, ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
