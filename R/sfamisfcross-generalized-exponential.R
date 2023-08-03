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
# Convolution: generalized genexponential - normal                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for misf genexpo-normal distribution
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
cmisfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A1 <- S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# cauchit specification class membership
cmisfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A1 <- S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# probit specification class membership
cmisfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A1 <- S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# cloglog specification class membership
cmisfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  A1 <- S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for misf genexpo-normal distribution
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
cstmisfgenexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- cstgenexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGenExpo <- NULL
  } else {
    cat("Initialization: SFA + generalized genexponential - normal distributions...\n")
    initGenExpo <- maxLik::maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradgenexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initGenExpo$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1),
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("MISF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Gradient of the likelihood function ----------
#' gradient for misf genexpo-normal distribution
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
cgradmisfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv <- exp(Wv)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv_h/ewu1_h + S * (epsilon)/ewv_h)
  musig2 <- (ewv_h/ewu2_h + S * (epsilon)/ewv_h)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  epsi1 <- (2 * (ewv_h/ewu1_h) + S * (epsilon)/ewv_h)
  epsi2 <- (2 * (ewv_h/ewu2_h) + S * (epsilon)/ewv_h)
  depsi1 <- dnorm(-epsi1, 0, 1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  pepsi2 <- pnorm(-epsi2)
  sigx1_1 <- (dmusig1/ewv_h - pmusig1/ewu1_h)
  sigx1_2 <- (dmusig2/ewv_h - pmusig2/ewu2_h)
  sigx2_1 <- (depsi1/ewv_h - 2 * (pepsi1/ewu1_h))
  sigx2_2 <- (depsi2/ewv_h - 2 * (pepsi2/ewu2_h))
  sigx3_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx3_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx4_1 <- exp(2 * (ewv/ewu1) + 2 * (S * (epsilon)/ewu1_h))
  sigx4_2 <- exp(2 * (ewv/ewu2) + 2 * (S * (epsilon)/ewu2_h))
  sigx5_1 <- (sigx1_1 * sigx3_1 - sigx2_1 * sigx4_1)
  sigx5_2 <- (sigx1_2 * sigx3_2 - sigx2_2 * sigx4_2)
  sigx6_1 <- (sigx5_1 * ewz/(wzdeno * ewu1_h))
  sigx6_2 <- (sigx5_2 * prC/ewu2_h)
  sigx7_1 <- (sigx3_1 * pmusig1 - sigx4_1 * pepsi1)
  sigx7_2 <- (sigx3_2 * pmusig2 - sigx4_2 * pepsi2)
  sigx8_1 <- (sigx7_1 * ewz/(wzdeno * ewu1_h))
  sigx8_2 <- (prC * sigx7_2/ewu2_h)
  s8 <- (2 * sigx8_2 + 2 * sigx8_1)
  sigx9 <- (2 * sigx6_1 + 2 * sigx6_2)/s8
  sigx10_1 <- (dmusig1 * ewv_h/ewu1_h)
  sigx10_2 <- (dmusig2 * ewv_h/ewu2_h)
  sigx11_1 <- (0.5 * (S * (epsilon)/ewu1_h) + 2 * (ewu1 * ewv/(2 * ewu1)^2))
  sigx11_2 <- (0.5 * (S * (epsilon)/ewu2_h) + 2 * (ewu2 * ewv/(2 * ewu2)^2))
  uvsi1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewu1_h)
  uvsi2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewu2_h)
  sigx12_1 <- (depsi1 * ewv_h/ewu1_h - uvsi1 * pepsi1)
  sigx12_2 <- (depsi2 * ewv_h/ewu2_h - uvsi2 * pepsi2)
  sp1 <- (0.5 * sigx10_1 - sigx11_1 * pmusig1)
  sp2 <- (0.5 * sigx10_2 - sigx11_2 * pmusig2)
  sigx13_1 <- (sp1 * sigx3_1 - sigx12_1 * sigx4_1)
  sigx13_2 <- (sp2 * sigx3_2 - (sigx12_2 * sigx4_2 + 0.5 * sigx7_2))
  sigx14_1 <- (wzdeno * ewu1_h * sigx7_1/(wzdeno * ewu1_h)^2)
  sigx15_1 <- (sigx13_1/(wzdeno * ewu1_h) - 0.5 * sigx14_1)
  sigx16_1 <- (0.5 * (ewv_h/ewu1_h) - 0.5 * (S * (epsilon)/ewv_h))
  sigx16_2 <- (0.5 * (ewv_h/ewu2_h) - 0.5 * (S * (epsilon)/ewv_h))
  sigx17_1 <- (ewv * pmusig1/(2 * ewu1) - sigx16_1 * dmusig1)
  sigx17_2 <- (ewv * pmusig2/(2 * ewu2) - sigx16_2 * dmusig2)
  sigx18_1 <- (ewv_h/ewu1_h - 0.5 * (S * (epsilon)/ewv_h))
  sigx18_2 <- (ewv_h/ewu2_h - 0.5 * (S * (epsilon)/ewv_h))
  sigx19_1 <- (sigx3_1 * sigx17_1 - (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx18_1) *
    sigx4_1)
  sigx19_2 <- (sigx3_2 * sigx17_2 - (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx18_2) *
    sigx4_2)
  sigx20_1 <- (1/(wzdeno * ewu1_h) - ewu1_h * ewz/(wzdeno * ewu1_h)^2)
  sigx20_2 <- (prC * sigx7_2/(wzdeno * ewu2_h))
  sigx21 <- (2 * (sigx20_1 * sigx7_1) - 2 * sigx20_2)
  wzsig <- (wzdeno * s8 * ewu1_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx9, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = 2 * (sigx15_1 * ewz/s8), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = 2 * (sigx13_2 * prC/(s8 * ewu2_h)), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = 2 * (sigx19_1 * ewz/wzsig) + 2 * (prC * sigx19_2/(s8 * ewu2_h)),
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx21 * ewz/s8, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cauchit specification class membership
cgradmisfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr <- exp(Wv/2)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr/ewusr1 + S * (epsilon)/ewvsr)
  musig2 <- (ewvsr/ewusr2 + S * (epsilon)/ewvsr)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewvsr/ewusr1) + S * (epsilon)/ewvsr)
  epsivu2 <- (2 * (ewvsr/ewusr2) + S * (epsilon)/ewvsr)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  dpuv1 <- (depsi1/ewvsr - 2 * (pepsi1/ewusr1))
  dpuv2 <- (depsi2/ewvsr - 2 * (pepsi2/ewusr2))
  euv1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewusr1)
  euv2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx0_1 <- exp(2 * (ewv/ewu1) + 2 * epsiu1)
  sigx0_2 <- exp(2 * (ewv/ewu2) + 2 * epsiu2)
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- ((sigx1_1 * sigx2_1 - dpuv1 * sigx0_1) * ewz1/ewusr1)
  sigx3_2 <- ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2) * ewz2/ewusr2)
  sigx14_1 <- (sigx2_1 * pmusig1 - sigx0_1 * pepsi1)
  sigx14_2 <- (sigx2_2 * pmusig2 - sigx0_2 * pepsi2)
  sigx4_1 <- (ewz1 * sigx14_1/ewusr1)
  sigx4_2 <- (ewz2 * sigx14_2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) * pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) * pmusig2)
  sigx8_1 <- (depsi1 * ewvsr/ewusr1 - euv1 * pepsi1)
  sigx8_2 <- (depsi2 * ewvsr/ewusr2 - euv2 * pepsi2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10_1 <- (sigx7_1 * sigx2_1 - (sigx8_1 * sigx0_1 + 0.5 * sigx14_1))
  sigx10_2 <- (sigx7_2 * sigx2_2 - (sigx8_2 * sigx0_2 + 0.5 * sigx14_2))
  sigx11_1 <- (ewvsr/ewusr1 - 0.5 * (S * (epsilon)/ewvsr))
  sigx11_2 <- (ewvsr/ewusr2 - 0.5 * (S * (epsilon)/ewvsr))
  sigx12_1 <- (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx11_1)
  sigx12_2 <- (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx11_2)
  sigx13_1 <- (sigx2_1 * sigx9_1 - sigx12_1 * sigx0_1)
  sigx13_2 <- (sigx2_2 * sigx9_2 - sigx12_2 * sigx0_2)
  sigx15 <- (pi * ((Wz)^2 + 1) * (2 * sigx4_2 + 2 * sigx4_1))
  sigx16 <- (2 * (ewz2 * sigx13_2/ewusr2) + 2 * (ewz1 * sigx13_1/ewusr1))
  sigx17 <- (2 * (sigx14_1/ewusr1) - 2 * (sigx14_2/ewusr2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (2 * sigx3_1 + 2 * sigx3_2)/(2 *
    sigx4_2 + 2 * sigx4_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx10_1 * ewz1/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1)), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = 2 * (sigx10_2 * ewz2/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2)),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16/(2 * sigx4_2 + 2 * sigx4_1),
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx17/sigx15, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# probit specification class membership
cgradmisfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr <- exp(Wv/2)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr/ewusr1 + S * (epsilon)/ewvsr)
  musig2 <- (ewvsr/ewusr2 + S * (epsilon)/ewvsr)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewvsr/ewusr1) + S * (epsilon)/ewvsr)
  epsivu2 <- (2 * (ewvsr/ewusr2) + S * (epsilon)/ewvsr)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  dpuv1 <- (depsi1/ewvsr - 2 * (pepsi1/ewusr1))
  dpuv2 <- (depsi2/ewvsr - 2 * (pepsi2/ewusr2))
  euv1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewusr1)
  euv2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx0_1 <- exp(2 * (ewv/ewu1) + 2 * epsiu1)
  sigx0_2 <- exp(2 * (ewv/ewu2) + 2 * epsiu2)
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- ((sigx1_1 * sigx2_1 - dpuv1 * sigx0_1) * pwZ/ewusr1)
  sigx3_2 <- ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2) * (1 - pwZ)/ewusr2)
  sigx14_1 <- (sigx2_1 * pmusig1 - sigx0_1 * pepsi1)
  sigx14_2 <- (sigx2_2 * pmusig2 - sigx0_2 * pepsi2)
  sigx4_1 <- (sigx14_1 * pwZ/ewusr1)
  sigx4_2 <- ((1 - pwZ) * sigx14_2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) * pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) * pmusig2)
  sigx8_1 <- (depsi1 * ewvsr/ewusr1 - euv1 * pepsi1)
  sigx8_2 <- (depsi2 * ewvsr/ewusr2 - euv2 * pepsi2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10_1 <- (sigx7_1 * sigx2_1 - (sigx8_1 * sigx0_1 + 0.5 * sigx14_1))
  sigx10_2 <- (sigx7_2 * sigx2_2 - (sigx8_2 * sigx0_2 + 0.5 * sigx14_2))
  sigx11_1 <- (ewvsr/ewusr1 - 0.5 * (S * (epsilon)/ewvsr))
  sigx11_2 <- (ewvsr/ewusr2 - 0.5 * (S * (epsilon)/ewvsr))
  sigx12_1 <- (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx11_1)
  sigx12_2 <- (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx11_2)
  sigx13_1 <- (sigx2_1 * sigx9_1 - sigx12_1 * sigx0_1)
  sigx13_2 <- (sigx2_2 * sigx9_2 - sigx12_2 * sigx0_2)
  sigx16 <- (2 * ((1 - pwZ) * sigx13_2/ewusr2) + 2 * (sigx13_1 * pwZ/ewusr1))
  sigx17 <- (2 * (sigx14_1/ewusr1) - 2 * (sigx14_2/ewusr2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (2 * sigx3_1 + 2 * sigx3_2)/(2 *
    sigx4_2 + 2 * sigx4_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx10_1 * pwZ/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1)), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = 2 * (sigx10_2 * (1 - pwZ)/((2 * sigx4_2 + 2 * sigx4_1) *
      ewusr2)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16/(2 * sigx4_2 +
    2 * sigx4_1), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx17 * dwZ/(2 *
    sigx4_2 + 2 * sigx4_1), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cloglog specification class membership
cgradmisfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr <- exp(Wv/2)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr/ewusr1 + S * (epsilon)/ewvsr)
  musig2 <- (ewvsr/ewusr2 + S * (epsilon)/ewvsr)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewvsr/ewusr1) + S * (epsilon)/ewvsr)
  epsivu2 <- (2 * (ewvsr/ewusr2) + S * (epsilon)/ewvsr)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  dpuv1 <- (depsi1/ewvsr - 2 * (pepsi1/ewusr1))
  dpuv2 <- (depsi2/ewvsr - 2 * (pepsi2/ewusr2))
  euv1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewusr1)
  euv2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx0_1 <- exp(2 * (ewv/ewu1) + 2 * epsiu1)
  sigx0_2 <- exp(2 * (ewv/ewu2) + 2 * epsiu2)
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- ((sigx1_1 * sigx2_1 - dpuv1 * sigx0_1) * (1 - prZ)/ewusr1)
  sigx3_2 <- ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2) * prZ/ewusr2)
  sigx14_1 <- (sigx2_1 * pmusig1 - sigx0_1 * pepsi1)
  sigx14_2 <- (sigx2_2 * pmusig2 - sigx0_2 * pepsi2)
  sigx4_1 <- ((1 - prZ) * sigx14_1/ewusr1)
  sigx4_2 <- (prZ * sigx14_2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) * pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) * pmusig2)
  sigx8_1 <- (depsi1 * ewvsr/ewusr1 - euv1 * pepsi1)
  sigx8_2 <- (depsi2 * ewvsr/ewusr2 - euv2 * pepsi2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10_1 <- (sigx7_1 * sigx2_1 - (sigx8_1 * sigx0_1 + 0.5 * sigx14_1))
  sigx10_2 <- (sigx7_2 * sigx2_2 - (sigx8_2 * sigx0_2 + 0.5 * sigx14_2))
  sigx11_1 <- (ewvsr/ewusr1 - 0.5 * (S * (epsilon)/ewvsr))
  sigx11_2 <- (ewvsr/ewusr2 - 0.5 * (S * (epsilon)/ewvsr))
  sigx12_1 <- (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx11_1)
  sigx12_2 <- (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx11_2)
  sigx13_1 <- (sigx2_1 * sigx9_1 - sigx12_1 * sigx0_1)
  sigx13_2 <- (sigx2_2 * sigx9_2 - sigx12_2 * sigx0_2)
  sigx15 <- (pi * ((Wz)^2 + 1) * (2 * sigx4_2 + 2 * sigx4_1))
  sigx16 <- (2 * ((1 - prZ) * sigx13_1/ewusr1) + 2 * (prZ * sigx13_2/ewusr2))
  sigx17 <- (2 * (sigx14_1/ewusr1) - 2 * (sigx14_2/ewusr2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (2 * sigx3_1 + 2 * sigx3_2)/(2 *
    sigx4_2 + 2 * sigx4_1), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx10_1 * (1 - prZ)/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx10_2 * prZ/((2 * sigx4_2 + 2 *
      sigx4_1) * ewusr2)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16/(2 *
      sigx4_2 + 2 * sigx4_1), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx17 *
      prZ * ewz/(2 * sigx4_1 + 2 * sigx4_2), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for misf genexpo-normal distribution
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
chessmisfgenexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewv <- exp(Wv)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv_h/ewu1_h + S * (epsilon)/ewv_h)
  musig2 <- (ewv_h/ewu2_h + S * (epsilon)/ewv_h)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  epsi1 <- (2 * (ewv_h/ewu1_h) + S * (epsilon)/ewv_h)
  epsi2 <- (2 * (ewv_h/ewu2_h) + S * (epsilon)/ewv_h)
  depsi1 <- dnorm(-epsi1, 0, 1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  pepsi2 <- pnorm(-epsi2)
  sigx1_1 <- (dmusig1/ewv_h - pmusig1/ewu1_h)
  sigx1_2 <- (dmusig2/ewv_h - pmusig2/ewu2_h)
  sigx2_1 <- (depsi1/ewv_h - 2 * (pepsi1/ewu1_h))
  sigx2_2 <- (depsi2/ewv_h - 2 * (pepsi2/ewu2_h))
  sigx3_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx3_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx4_1 <- exp(2 * (ewv/ewu1) + 2 * (S * (epsilon)/ewu1_h))
  sigx4_2 <- exp(2 * (ewv/ewu2) + 2 * (S * (epsilon)/ewu2_h))
  sigx5_1 <- (sigx1_1 * sigx3_1 - sigx2_1 * sigx4_1)
  sigx5_2 <- (sigx1_2 * sigx3_2 - sigx2_2 * sigx4_2)
  sigx6_1 <- (sigx5_1 * ewz/(wzdeno * ewu1_h))
  sigx6_2 <- (sigx5_2 * prC/ewu2_h)
  sigx7_1 <- (sigx3_1 * pmusig1 - sigx4_1 * pepsi1)
  sigx7_2 <- (sigx3_2 * pmusig2 - sigx4_2 * pepsi2)
  sigx8_1 <- (sigx7_1 * ewz/(wzdeno * ewu1_h))
  sigx8_2 <- (prC * sigx7_2/ewu2_h)
  s8 <- (2 * sigx8_2 + 2 * sigx8_1)
  sigx9 <- (2 * sigx6_1 + 2 * sigx6_2)/s8
  sigx10_1 <- (dmusig1 * ewv_h/ewu1_h)
  sigx10_2 <- (dmusig2 * ewv_h/ewu2_h)
  sigx11_1 <- (0.5 * (S * (epsilon)/ewu1_h) + 2 * (ewu1 * ewv/(2 * ewu1)^2))
  sigx11_2 <- (0.5 * (S * (epsilon)/ewu2_h) + 2 * (ewu2 * ewv/(2 * ewu2)^2))
  uvsi1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewu1_h)
  uvsi2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewu2_h)
  sigx12_1 <- (depsi1 * ewv_h/ewu1_h - uvsi1 * pepsi1)
  sigx12_2 <- (depsi2 * ewv_h/ewu2_h - uvsi2 * pepsi2)
  sp1 <- (0.5 * sigx10_1 - sigx11_1 * pmusig1)
  sp2 <- (0.5 * sigx10_2 - sigx11_2 * pmusig2)
  sigx13_1 <- (sp1 * sigx3_1 - sigx12_1 * sigx4_1)
  sigx13_2 <- (sp2 * sigx3_2 - (sigx12_2 * sigx4_2 + 0.5 * sigx7_2))
  sigx14_1 <- (wzdeno * ewu1_h * sigx7_1/(wzdeno * ewu1_h)^2)
  sigx15_1 <- (sigx13_1/(wzdeno * ewu1_h) - 0.5 * sigx14_1)
  sigx16_1 <- (0.5 * (ewv_h/ewu1_h) - 0.5 * (S * (epsilon)/ewv_h))
  sigx16_2 <- (0.5 * (ewv_h/ewu2_h) - 0.5 * (S * (epsilon)/ewv_h))
  sigx17_1 <- (ewv * pmusig1/(2 * ewu1) - sigx16_1 * dmusig1)
  sigx17_2 <- (ewv * pmusig2/(2 * ewu2) - sigx16_2 * dmusig2)
  sigx18_1 <- (ewv_h/ewu1_h - 0.5 * (S * (epsilon)/ewv_h))
  sigx18_2 <- (ewv_h/ewu2_h - 0.5 * (S * (epsilon)/ewv_h))
  sigx19_1 <- (sigx3_1 * sigx17_1 - (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx18_1) *
    sigx4_1)
  sigx19_2 <- (sigx3_2 * sigx17_2 - (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx18_2) *
    sigx4_2)
  sigx20_1 <- (1/(wzdeno * ewu1_h) - ewu1_h * ewz/(wzdeno * ewu1_h)^2)
  sigx20_2 <- (prC * sigx7_2/(wzdeno * ewu2_h))
  sigx21 <- (2 * (sigx20_1 * sigx7_1) - 2 * sigx20_2)
  sigx22_1 <- (0.5 - wzdeno^2 * ewu1_h^2/(wzdeno * ewu1_h)^2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * (2 * ((((musig1/ewv_h - 1/ewu1_h) * dmusig1/ewv_h - sigx1_1/ewu1_h) *
    sigx3_1 - ((epsi1/ewv_h - 2/ewu1_h) * depsi1/ewv_h - 2 * (sigx2_1/ewu1_h)) *
    sigx4_1) * ewz/(wzdeno * ewu1_h)) + 2 * ((((musig2/ewv_h - 1/ewu2_h) * dmusig2/ewv_h -
    sigx1_2/ewu2_h) * sigx3_2 - ((epsi2/ewv_h - 2/ewu2_h) * depsi2/ewv_h - 2 *
    (sigx2_2/ewu2_h)) * sigx4_2) * prC/ewu2_h) - (2 * sigx6_1 + 2 * sigx6_2)^2/s8)/s8,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewu1_h) + 2 * (ewu1 *
      ewv/(2 * ewu1)^2)) * pmusig1 - 0.5 * sigx10_1)/ewu1_h + (0.5 * (musig1/ewu1_h) -
      sigx11_1/ewv_h) * dmusig1) * sigx3_1 - ((epsi1/ewu1_h - uvsi1/ewv_h) *
      depsi1 + (pepsi1 - 2 * sigx12_1)/ewu1_h) * sigx4_1)/(wzdeno * ewu1_h) -
      (sigx15_1 * sigx9 + 0.5 * (sigx5_1 * wzdeno * ewu1_h/(wzdeno * ewu1_h)^2))) *
      ewz/s8), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * (S * (epsilon)/ewu2_h) +
      2 * (ewu2 * ewv/(2 * ewu2)^2)) * pmusig2 - 0.5 * sigx10_2)/ewu2_h + (0.5 *
      (musig2/ewu2_h) - sigx11_2/ewv_h) * dmusig2) * sigx3_2 - (((epsi2/ewu2_h -
      uvsi2/ewv_h) * depsi2 + (pepsi2 - 2 * sigx12_2)/ewu2_h) * sigx4_2 + 0.5 *
      sigx5_2))/(s8 * ewu2_h) - sigx13_2 * (2 * sigx6_1 + 2 * sigx6_2) * ewu2_h/(s8 *
      ewu2_h)^2) * prC), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((dmusig1 * (ewv/(2 * ewu1) - (sigx16_1 *
      musig1 + 0.5))/ewv_h - sigx17_1/ewu1_h) * sigx3_1 - ((2 * (ewv/ewu1) -
      (epsi1 * sigx18_1 + 0.5)) * depsi1/ewv_h - 2 * ((2 * (ewv * pepsi1/ewu1) -
      depsi1 * sigx18_1)/ewu1_h)) * sigx4_1) * ewz/(wzdeno * ewu1_h)) + 2 *
      (((dmusig2 * (ewv/(2 * ewu2) - (sigx16_2 * musig2 + 0.5))/ewv_h - sigx17_2/ewu2_h) *
        sigx3_2 - ((2 * (ewv/ewu2) - (epsi2 * sigx18_2 + 0.5)) * depsi2/ewv_h -
        2 * ((2 * (ewv * pepsi2/ewu2) - depsi2 * sigx18_2)/ewu2_h)) * sigx4_2) *
        prC/ewu2_h) - (2 * sigx6_1 + 2 * sigx6_2) * (2 * (prC * sigx19_2/ewu2_h) +
      2 * (sigx19_1 * ewz/(wzdeno * ewu1_h)))/s8)/s8, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * (sigx5_1 *
    sigx20_1) - ((2 * sigx6_1 + 2 * sigx6_2) * sigx21/s8 + 2 * (sigx5_2 * prC/(wzdeno *
    ewu2_h)))) * ewz/s8, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewv_h * musig1/ewu1_h) -
      0.5) - 0.5 * sigx11_1) * dmusig1 * ewv_h/ewu1_h - (sp1 * sigx11_1 + (2 *
      ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv/(2 * ewu1)^2) - 0.25 *
      (S * (epsilon)/ewu1_h)) * pmusig1)) * sigx3_1 - (((epsi1 * ewv_h - S *
      (epsilon))/ewu1_h - (0.5 + 2 * (ewv/ewu1))) * depsi1 * ewv_h/ewu1_h +
      (0.5 * (S * (epsilon)/ewu1_h) + 2 * (ewv/ewu1)) * pepsi1 - uvsi1 * sigx12_1) *
      sigx4_1)/(wzdeno * ewu1_h) - ((0.5 * (sigx22_1 * sigx7_1 + sp1 * sigx3_1 -
      sigx12_1 * sigx4_1) + 0.5 * sigx13_1) * wzdeno * ewu1_h/(wzdeno * ewu1_h)^2 +
      2 * (sigx15_1^2 * ewz/s8))) * ewz/s8), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (sigx15_1 * sigx13_2 * prC * ewu2_h * ewz/(s8 *
      ewu2_h)^2)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((((dmusig1 * ewv_h/(4 * (ewu1 * ewu1_h)) - 2 * (ewu1 * pmusig1/(2 *
      ewu1)^2)) * ewv - ((0.5 * (sigx16_1 * musig1) - 0.25) * dmusig1 * ewv_h/ewu1_h +
      sigx11_1 * sigx17_1)) * sigx3_1 - (2 * ((depsi1 * ewv_h/ewu1_h - pepsi1) *
      ewv/ewu1) - ((epsi1 * sigx18_1 - 0.5) * depsi1 * ewv_h/ewu1_h + (2 *
      (ewv * pepsi1/ewu1) - depsi1 * sigx18_1) * uvsi1)) * sigx4_1)/(wzdeno *
      ewu1_h) - 0.5 * (wzdeno * ewu1_h * sigx19_1/(wzdeno * ewu1_h)^2)) - 2 *
      (sigx15_1 * (2 * (prC * sigx19_2/ewu2_h) + 2 * (sigx19_1 * ewz/(wzdeno *
        ewu1_h)))/s8)) * ewz/s8, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * (sigx13_1 * sigx20_1 - (sigx22_1 * ewz + 0.5 * wzdeno) * ewu1_h * sigx7_1/(wzdeno *
      ewu1_h)^2) - 2 * (sigx15_1 * sigx21 * ewz/s8)) * ewz/s8, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((0.5 *
    (0.5 * (ewv_h * musig2/ewu2_h) - 0.5) - 0.5 * sigx11_2) * dmusig2 * ewv_h/ewu2_h -
    (sp2 * sigx11_2 + (2 * ((1 - 8 * (ewu2^2/(2 * ewu2)^2)) * ewu2 * ewv/(2 *
      ewu2)^2) - 0.25 * (S * (epsilon)/ewu2_h)) * pmusig2)) * sigx3_2 - ((((epsi2 *
    ewv_h - S * (epsilon))/ewu2_h - (0.5 + 2 * (ewv/ewu2))) * depsi2 * ewv_h/ewu2_h +
    (0.5 * (S * (epsilon)/ewu2_h) + 2 * (ewv/ewu2)) * pepsi2 - uvsi2 * sigx12_2) *
    sigx4_2 + 0.5 * (sp2 * sigx3_2 - sigx12_2 * sigx4_2)))/(s8 * ewu2_h) - sigx13_2 *
    (0.5 * (s8 * ewu2_h) + 2 * (sigx13_2 * prC))/(s8 * ewu2_h)^2) * prC), FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    prC * (2 * (((dmusig2 * ewv_h/(4 * (ewu2 * ewu2_h)) - 2 * (ewu2 * pmusig2/(2 *
    ewu2)^2)) * ewv - ((0.5 * (sigx16_2 * musig2) - 0.25) * dmusig2 * ewv_h/ewu2_h +
    sigx11_2 * sigx17_2)) * sigx3_2 - ((2 * ((depsi2 * ewv_h/ewu2_h - pepsi2) *
    ewv/ewu2) - ((epsi2 * sigx18_2 - 0.5) * depsi2 * ewv_h/ewu2_h + (2 * (ewv *
    pepsi2/ewu2) - depsi2 * sigx18_2) * uvsi2)) * sigx4_2 + 0.5 * sigx19_2)) -
    2 * (sigx13_2 * (2 * (prC * sigx19_2/ewu2_h) + 2 * (sigx19_1 * ewz/(wzdeno *
      ewu1_h)))/s8))/(s8 * ewu2_h), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (prC * (2 * (sigx13_2 * sigx21/(s8 * ewu2_h)) + 2 * ((sp2 *
      sigx3_2 - sigx12_2 * sigx4_2)/(wzdeno * ewu2_h) - 0.5 * (wzdeno * ewu2_h *
      sigx7_2/(wzdeno * ewu2_h)^2))) * ewz/s8), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (2 * ((((sigx17_1/2 + (pmusig1 - sigx16_1 * dmusig1)/2) *
      ewv/ewu1 - (0.25 * (ewv_h/ewu1_h) + 0.25 * (S * (epsilon)/ewv_h) - sigx16_1^2 *
      musig1) * dmusig1) * sigx3_1 - ((2 * (2 * (ewv * pepsi1/ewu1) - depsi1 *
      sigx18_1) + 2 * (pepsi1 - depsi1 * sigx18_1)) * ewv/ewu1 - (0.25 * (S *
      (epsilon)/ewv_h) + 0.5 * (ewv_h/ewu1_h) - epsi1 * sigx18_1^2) * depsi1) *
      sigx4_1) * ewz/(wzdeno * ewu1_h)) + 2 * ((((sigx17_2/2 + (pmusig2 - sigx16_2 *
      dmusig2)/2) * ewv/ewu2 - (0.25 * (ewv_h/ewu2_h) + 0.25 * (S * (epsilon)/ewv_h) -
      sigx16_2^2 * musig2) * dmusig2) * sigx3_2 - ((2 * (2 * (ewv * pepsi2/ewu2) -
      depsi2 * sigx18_2) + 2 * (pepsi2 - depsi2 * sigx18_2)) * ewv/ewu2 - (0.25 *
      (S * (epsilon)/ewv_h) + 0.5 * (ewv_h/ewu2_h) - epsi2 * sigx18_2^2) *
      depsi2) * sigx4_2) * prC/ewu2_h) - (2 * (prC * sigx19_2/ewu2_h) + 2 *
      (sigx19_1 * ewz/(wzdeno * ewu1_h)))^2/s8)/s8, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx20_1 * sigx19_1) - ((2 * (prC * sigx19_2/ewu2_h) +
      2 * (sigx19_1 * ewz/(wzdeno * ewu1_h))) * sigx21/s8 + 2 * (prC * sigx19_2/(wzdeno *
      ewu2_h)))) * ewz/s8, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((2 * (prC * (1/(wzdeno^2 * ewu2_h) + ewu2_h/(wzdeno *
      ewu2_h)^2) * sigx7_2) - (sigx21^2/s8 + 2 * ((2 - 2 * (wzdeno * ewu1_h^2 *
      ewz/(wzdeno * ewu1_h)^2)) * ewu1_h * sigx7_1/(wzdeno * ewu1_h)^2))) *
      ewz + 2 * (sigx20_1 * sigx7_1) - 2 * sigx20_2) * ewz/s8, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cauchit specification class membership
chessmisfgenexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr <- exp(Wv/2)
  ewv <- exp(Wv)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr/ewusr1 + S * (epsilon)/ewvsr)
  musig2 <- (ewvsr/ewusr2 + S * (epsilon)/ewvsr)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewvsr/ewusr1) + S * (epsilon)/ewvsr)
  epsivu2 <- (2 * (ewvsr/ewusr2) + S * (epsilon)/ewvsr)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  dpuv1 <- (depsi1/ewvsr - 2 * (pepsi1/ewusr1))
  dpuv2 <- (depsi2/ewvsr - 2 * (pepsi2/ewusr2))
  euv1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewusr1)
  euv2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx0_1 <- exp(2 * (ewv/ewu1) + 2 * epsiu1)
  sigx0_2 <- exp(2 * (ewv/ewu2) + 2 * epsiu2)
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- ((sigx1_1 * sigx2_1 - dpuv1 * sigx0_1) * ewz1/ewusr1)
  sigx3_2 <- ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2) * ewz2/ewusr2)
  sigx14_1 <- (sigx2_1 * pmusig1 - sigx0_1 * pepsi1)
  sigx14_2 <- (sigx2_2 * pmusig2 - sigx0_2 * pepsi2)
  sigx4_1 <- (ewz1 * sigx14_1/ewusr1)
  sigx4_2 <- (ewz2 * sigx14_2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) * pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) * pmusig2)
  sigx8_1 <- (depsi1 * ewvsr/ewusr1 - euv1 * pepsi1)
  sigx8_2 <- (depsi2 * ewvsr/ewusr2 - euv2 * pepsi2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10_1 <- (sigx7_1 * sigx2_1 - (sigx8_1 * sigx0_1 + 0.5 * sigx14_1))
  sigx10_2 <- (sigx7_2 * sigx2_2 - (sigx8_2 * sigx0_2 + 0.5 * sigx14_2))
  sigx11_1 <- (ewvsr/ewusr1 - 0.5 * (S * (epsilon)/ewvsr))
  sigx11_2 <- (ewvsr/ewusr2 - 0.5 * (S * (epsilon)/ewvsr))
  sigx12_1 <- (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx11_1)
  sigx12_2 <- (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx11_2)
  sigx13_1 <- (sigx2_1 * sigx9_1 - sigx12_1 * sigx0_1)
  sigx13_2 <- (sigx2_2 * sigx9_2 - sigx12_2 * sigx0_2)
  sigx15 <- (pi * ((Wz)^2 + 1) * (2 * sigx4_2 + 2 * sigx4_1))
  sigx16 <- (2 * (ewz2 * sigx13_2/ewusr2) + 2 * (ewz1 * sigx13_1/ewusr1))
  sigx17 <- (2 * (sigx14_1/ewusr1) - 2 * (sigx14_2/ewusr2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (2 * ((((musig1/ewvsr - 1/ewusr1) * dmusig1/ewvsr - sigx1_1/ewusr1) * sigx2_1 -
      ((epsivu1/ewvsr - 2/ewusr1) * depsi1/ewvsr - 2 * (dpuv1/ewusr1)) * sigx0_1) *
      ewz1/ewusr1) + 2 * ((((musig2/ewvsr - 1/ewusr2) * dmusig2/ewvsr - sigx1_2/ewusr2) *
      sigx2_2 - ((epsivu2/ewvsr - 2/ewusr2) * depsi2/ewvsr - 2 * (dpuv2/ewusr2)) *
      sigx0_2) * ewz2/ewusr2) - (2 * sigx3_1 + 2 * sigx3_2)^2/(2 * sigx4_2 +
      2 * sigx4_1))/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * epsiu1 + 2 * sigx6_1) * pmusig1 -
      0.5 * sigx5_1)/ewusr1 + (0.5 * (musig1/ewusr1) - (0.5 * epsiu1 + 2 *
      sigx6_1)/ewvsr) * dmusig1) * sigx2_1 - (((epsivu1/ewusr1 - euv1/ewvsr) *
      depsi1 + (pepsi1 - 2 * sigx8_1)/ewusr1) * sigx0_1 + 0.5 * (sigx1_1 *
      sigx2_1 - dpuv1 * sigx0_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1) -
      sigx10_1 * (2 * sigx3_1 + 2 * sigx3_2) * ewusr1/((2 * sigx4_2 + 2 * sigx4_1) *
        ewusr1)^2) * ewz1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
      pmusig2 - 0.5 * sigx5_2)/ewusr2 + (0.5 * (musig2/ewusr2) - (0.5 * epsiu2 +
      2 * sigx6_2)/ewvsr) * dmusig2) * sigx2_2 - (((epsivu2/ewusr2 - euv2/ewvsr) *
      depsi2 + (pepsi2 - 2 * sigx8_2)/ewusr2) * sigx0_2 + 0.5 * (sigx1_2 *
      sigx2_2 - dpuv2 * sigx0_2)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2) -
      sigx10_2 * (2 * sigx3_1 + 2 * sigx3_2) * ewusr2/((2 * sigx4_2 + 2 * sigx4_1) *
        ewusr2)^2) * ewz2), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((dmusig1 * (ewv/(2 * ewu1) - (epsiv1 *
      musig1 + 0.5))/ewvsr - sigx9_1/ewusr1) * sigx2_1 - ((2 * (ewv/ewu1) -
      (epsivu1 * sigx11_1 + 0.5)) * depsi1/ewvsr - 2 * (sigx12_1/ewusr1)) *
      sigx0_1) * ewz1/ewusr1) + 2 * (((dmusig2 * (ewv/(2 * ewu2) - (epsiv2 *
      musig2 + 0.5))/ewvsr - sigx9_2/ewusr2) * sigx2_2 - ((2 * (ewv/ewu2) -
      (epsivu2 * sigx11_2 + 0.5)) * depsi2/ewvsr - 2 * (sigx12_2/ewusr2)) *
      sigx0_2) * ewz2/ewusr2) - (2 * sigx3_1 + 2 * sigx3_2) * sigx16/(2 * sigx4_2 +
      2 * sigx4_1))/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((2 * ((sigx1_1 *
    sigx2_1 - dpuv1 * sigx0_1)/ewusr1) - 2 * ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2)/ewusr2))/sigx15 -
    pi * ((Wz)^2 + 1) * (2 * sigx3_1 + 2 * sigx3_2) * sigx17/sigx15^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewvsr * musig1/ewusr1) -
      0.5) - 0.5 * (0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr/ewusr1 -
      (sigx7_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 - 8 * (ewu1^2/(2 *
        ewu1)^2)) * ewu1 * ewv/(2 * ewu1)^2) - 0.25 * epsiu1) * pmusig1)) *
      sigx2_1 - ((((epsivu1 * ewvsr - S * (epsilon))/ewusr1 - (0.5 + 2 * (ewv/ewu1))) *
      depsi1 * ewvsr/ewusr1 + (0.5 * epsiu1 + 2 * (ewv/ewu1)) * pepsi1 - euv1 *
      sigx8_1) * sigx0_1 + 0.5 * (sigx7_1 * sigx2_1 - sigx8_1 * sigx0_1)))/((2 *
      sigx4_2 + 2 * sigx4_1) * ewusr1) - sigx10_1 * (0.5 * ((2 * sigx4_2 +
      2 * sigx4_1) * ewusr1) + 2 * (sigx10_1 * ewz1))/((2 * sigx4_2 + 2 * sigx4_1) *
      ewusr1)^2) * ewz1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (sigx10_1 * sigx10_2 * ewz2 * ewz1 * ewusr2/(((2 *
      sigx4_2 + 2 * sigx4_1) * ewusr2)^2 * ewusr1))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ewz1 * (2 * (((dmusig1 * ewvsr/(4 * (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 *
    ewu1)^2)) * ewv - ((0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 * ewvsr/ewusr1 +
    (0.5 * epsiu1 + 2 * sigx6_1) * sigx9_1)) * sigx2_1 - ((2 * ((depsi1 * ewvsr/ewusr1 -
    pepsi1) * ewv/ewu1) - ((epsivu1 * sigx11_1 - 0.5) * depsi1 * ewvsr/ewusr1 +
    sigx12_1 * euv1)) * sigx0_1 + 0.5 * sigx13_1)) - 2 * (sigx10_1 * sigx16/(2 *
    sigx4_2 + 2 * sigx4_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx10_1 * (2/sigx15 - 2 * (pi * ((Wz)^2 + 1) * ewz1 * sigx17/sigx15^2))/ewusr1,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((0.5 *
    (0.5 * (ewvsr * musig2/ewusr2) - 0.5) - 0.5 * (0.5 * epsiu2 + 2 * sigx6_2)) *
    dmusig2 * ewvsr/ewusr2 - (sigx7_2 * (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 -
    8 * (ewu2^2/(2 * ewu2)^2)) * ewu2 * ewv/(2 * ewu2)^2) - 0.25 * epsiu2) *
    pmusig2)) * sigx2_2 - ((((epsivu2 * ewvsr - S * (epsilon))/ewusr2 - (0.5 +
    2 * (ewv/ewu2))) * depsi2 * ewvsr/ewusr2 + (0.5 * epsiu2 + 2 * (ewv/ewu2)) *
    pepsi2 - euv2 * sigx8_2) * sigx0_2 + 0.5 * (sigx7_2 * sigx2_2 - sigx8_2 *
    sigx0_2)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2) - sigx10_2 * (0.5 * ((2 *
    sigx4_2 + 2 * sigx4_1) * ewusr2) + 2 * (sigx10_2 * ewz2))/((2 * sigx4_2 +
    2 * sigx4_1) * ewusr2)^2) * ewz2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ewz2 * (2 * (((dmusig2 * ewvsr/(4 * (ewu2 * ewusr2)) - 2 * (ewu2 * pmusig2/(2 *
    ewu2)^2)) * ewv - ((0.5 * (epsiv2 * musig2) - 0.25) * dmusig2 * ewvsr/ewusr2 +
    (0.5 * epsiu2 + 2 * sigx6_2) * sigx9_2)) * sigx2_2 - ((2 * ((depsi2 * ewvsr/ewusr2 -
    pepsi2) * ewv/ewu2) - ((epsivu2 * sigx11_2 - 0.5) * depsi2 * ewvsr/ewusr2 +
    sigx12_2 * euv2)) * sigx0_2 + 0.5 * sigx13_2)) - 2 * (sigx10_2 * sigx16/(2 *
    sigx4_2 + 2 * sigx4_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx10_2 * (2 * (pi * ((Wz)^2 + 1) * ewz2 * sigx17/sigx15^2) +
      2/sigx15)/ewusr2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (2 * ((((sigx9_1/2 + (pmusig1 - epsiv1 * dmusig1)/2) * ewv/ewu1 -
      (0.25 * (ewvsr/ewusr1) + 0.25 * (S * (epsilon)/ewvsr) - epsiv1^2 * musig1) *
        dmusig1) * sigx2_1 - ((2 * sigx12_1 + 2 * (pepsi1 - depsi1 * sigx11_1)) *
      ewv/ewu1 - (0.25 * (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr1) - epsivu1 *
      sigx11_1^2) * depsi1) * sigx0_1) * ewz1/ewusr1) + 2 * ((((sigx9_2/2 +
      (pmusig2 - epsiv2 * dmusig2)/2) * ewv/ewu2 - (0.25 * (ewvsr/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr) - epsiv2^2 * musig2) * dmusig2) * sigx2_2 -
      ((2 * sigx12_2 + 2 * (pepsi2 - depsi2 * sigx11_2)) * ewv/ewu2 - (0.25 *
        (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr2) - epsivu2 * sigx11_2^2) *
        depsi2) * sigx0_2) * ewz2/ewusr2) - sigx16^2/(2 * sigx4_2 + 2 * sigx4_1))/(2 *
      sigx4_2 + 2 * sigx4_1), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx13_1/ewusr1) - 2 * (sigx13_2/ewusr2))/sigx15 -
      pi * ((Wz)^2 + 1) * sigx16 * sigx17/sigx15^2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (sigx17 * (2 * (sigx14_1/ewusr1) + 2 * (pi *
      Wz * (2 * sigx4_2 + 2 * sigx4_1)) - 2 * (sigx14_2/ewusr2))/sigx15^2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# probit specification class membership
chessmisfgenexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr <- exp(Wv/2)
  ewv <- exp(Wv)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr/ewusr1 + S * (epsilon)/ewvsr)
  musig2 <- (ewvsr/ewusr2 + S * (epsilon)/ewvsr)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewvsr/ewusr1) + S * (epsilon)/ewvsr)
  epsivu2 <- (2 * (ewvsr/ewusr2) + S * (epsilon)/ewvsr)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  dpuv1 <- (depsi1/ewvsr - 2 * (pepsi1/ewusr1))
  dpuv2 <- (depsi2/ewvsr - 2 * (pepsi2/ewusr2))
  euv1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewusr1)
  euv2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx0_1 <- exp(2 * (ewv/ewu1) + 2 * epsiu1)
  sigx0_2 <- exp(2 * (ewv/ewu2) + 2 * epsiu2)
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- ((sigx1_1 * sigx2_1 - dpuv1 * sigx0_1) * pwZ/ewusr1)
  sigx3_2 <- ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2) * (1 - pwZ)/ewusr2)
  sigx14_1 <- (sigx2_1 * pmusig1 - sigx0_1 * pepsi1)
  sigx14_2 <- (sigx2_2 * pmusig2 - sigx0_2 * pepsi2)
  sigx4_1 <- (sigx14_1 * pwZ/ewusr1)
  sigx4_2 <- ((1 - pwZ) * sigx14_2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) * pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) * pmusig2)
  sigx8_1 <- (depsi1 * ewvsr/ewusr1 - euv1 * pepsi1)
  sigx8_2 <- (depsi2 * ewvsr/ewusr2 - euv2 * pepsi2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10_1 <- (sigx7_1 * sigx2_1 - (sigx8_1 * sigx0_1 + 0.5 * sigx14_1))
  sigx10_2 <- (sigx7_2 * sigx2_2 - (sigx8_2 * sigx0_2 + 0.5 * sigx14_2))
  sigx11_1 <- (ewvsr/ewusr1 - 0.5 * (S * (epsilon)/ewvsr))
  sigx11_2 <- (ewvsr/ewusr2 - 0.5 * (S * (epsilon)/ewvsr))
  sigx12_1 <- (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx11_1)
  sigx12_2 <- (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx11_2)
  sigx13_1 <- (sigx2_1 * sigx9_1 - sigx12_1 * sigx0_1)
  sigx13_2 <- (sigx2_2 * sigx9_2 - sigx12_2 * sigx0_2)
  sigx16 <- (2 * ((1 - pwZ) * sigx13_2/ewusr2) + 2 * (sigx13_1 * pwZ/ewusr1))
  sigx17 <- (2 * (sigx14_1/ewusr1) - 2 * (sigx14_2/ewusr2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (2 * ((((musig1/ewvsr - 1/ewusr1) * dmusig1/ewvsr - sigx1_1/ewusr1) * sigx2_1 -
      ((epsivu1/ewvsr - 2/ewusr1) * depsi1/ewvsr - 2 * (dpuv1/ewusr1)) * sigx0_1) *
      pwZ/ewusr1) + 2 * ((((musig2/ewvsr - 1/ewusr2) * dmusig2/ewvsr - sigx1_2/ewusr2) *
      sigx2_2 - ((epsivu2/ewvsr - 2/ewusr2) * depsi2/ewvsr - 2 * (dpuv2/ewusr2)) *
      sigx0_2) * (1 - pwZ)/ewusr2) - (2 * sigx3_1 + 2 * sigx3_2)^2/(2 * sigx4_2 +
      2 * sigx4_1))/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * epsiu1 + 2 * sigx6_1) * pmusig1 -
      0.5 * sigx5_1)/ewusr1 + (0.5 * (musig1/ewusr1) - (0.5 * epsiu1 + 2 *
      sigx6_1)/ewvsr) * dmusig1) * sigx2_1 - (((epsivu1/ewusr1 - euv1/ewvsr) *
      depsi1 + (pepsi1 - 2 * sigx8_1)/ewusr1) * sigx0_1 + 0.5 * (sigx1_1 *
      sigx2_1 - dpuv1 * sigx0_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1) -
      sigx10_1 * (2 * sigx3_1 + 2 * sigx3_2) * ewusr1/((2 * sigx4_2 + 2 * sigx4_1) *
        ewusr1)^2) * pwZ), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
      pmusig2 - 0.5 * sigx5_2)/ewusr2 + (0.5 * (musig2/ewusr2) - (0.5 * epsiu2 +
      2 * sigx6_2)/ewvsr) * dmusig2) * sigx2_2 - (((epsivu2/ewusr2 - euv2/ewvsr) *
      depsi2 + (pepsi2 - 2 * sigx8_2)/ewusr2) * sigx0_2 + 0.5 * (sigx1_2 *
      sigx2_2 - dpuv2 * sigx0_2)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2) -
      sigx10_2 * (2 * sigx3_1 + 2 * sigx3_2) * ewusr2/((2 * sigx4_2 + 2 * sigx4_1) *
        ewusr2)^2) * (1 - pwZ)), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((dmusig1 * (ewv/(2 * ewu1) - (epsiv1 *
      musig1 + 0.5))/ewvsr - sigx9_1/ewusr1) * sigx2_1 - ((2 * (ewv/ewu1) -
      (epsivu1 * sigx11_1 + 0.5)) * depsi1/ewvsr - 2 * (sigx12_1/ewusr1)) *
      sigx0_1) * pwZ/ewusr1) + 2 * (((dmusig2 * (ewv/(2 * ewu2) - (epsiv2 *
      musig2 + 0.5))/ewvsr - sigx9_2/ewusr2) * sigx2_2 - ((2 * (ewv/ewu2) -
      (epsivu2 * sigx11_2 + 0.5)) * depsi2/ewvsr - 2 * (sigx12_2/ewusr2)) *
      sigx0_2) * (1 - pwZ)/ewusr2) - (2 * sigx3_1 + 2 * sigx3_2) * sigx16/(2 *
      sigx4_2 + 2 * sigx4_1))/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * ((sigx1_1 *
    sigx2_1 - dpuv1 * sigx0_1)/ewusr1) - ((2 * sigx3_1 + 2 * sigx3_2) * sigx17/(2 *
    sigx4_2 + 2 * sigx4_1) + 2 * ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2)/ewusr2))) *
    dwZ/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewvsr * musig1/ewusr1) -
      0.5) - 0.5 * (0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr/ewusr1 -
      (sigx7_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 - 8 * (ewu1^2/(2 *
        ewu1)^2)) * ewu1 * ewv/(2 * ewu1)^2) - 0.25 * epsiu1) * pmusig1)) *
      sigx2_1 - ((((epsivu1 * ewvsr - S * (epsilon))/ewusr1 - (0.5 + 2 * (ewv/ewu1))) *
      depsi1 * ewvsr/ewusr1 + (0.5 * epsiu1 + 2 * (ewv/ewu1)) * pepsi1 - euv1 *
      sigx8_1) * sigx0_1 + 0.5 * (sigx7_1 * sigx2_1 - sigx8_1 * sigx0_1)))/((2 *
      sigx4_2 + 2 * sigx4_1) * ewusr1) - sigx10_1 * (0.5 * ((2 * sigx4_2 +
      2 * sigx4_1) * ewusr1) + 2 * (sigx10_1 * pwZ))/((2 * sigx4_2 + 2 * sigx4_1) *
      ewusr1)^2) * pwZ), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (sigx10_1 * sigx10_2 * (1 - pwZ) * pwZ *
      ewusr2/(((2 * sigx4_2 + 2 * sigx4_1) * ewusr2)^2 * ewusr1))), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    pwZ * (2 * (((dmusig1 * ewvsr/(4 * (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 *
    ewu1)^2)) * ewv - ((0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 * ewvsr/ewusr1 +
    (0.5 * epsiu1 + 2 * sigx6_1) * sigx9_1)) * sigx2_1 - ((2 * ((depsi1 * ewvsr/ewusr1 -
    pepsi1) * ewv/ewu1) - ((epsivu1 * sigx11_1 - 0.5) * depsi1 * ewvsr/ewusr1 +
    sigx12_1 * euv1)) * sigx0_1 + 0.5 * sigx13_1)) - 2 * (sigx10_1 * sigx16/(2 *
    sigx4_2 + 2 * sigx4_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx10_1 * (2 - 2 * (sigx17 * pwZ/(2 * sigx4_2 + 2 * sigx4_1))) * dwZ/((2 *
    sigx4_2 + 2 * sigx4_1) * ewusr1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((0.5 *
    (0.5 * (ewvsr * musig2/ewusr2) - 0.5) - 0.5 * (0.5 * epsiu2 + 2 * sigx6_2)) *
    dmusig2 * ewvsr/ewusr2 - (sigx7_2 * (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 -
    8 * (ewu2^2/(2 * ewu2)^2)) * ewu2 * ewv/(2 * ewu2)^2) - 0.25 * epsiu2) *
    pmusig2)) * sigx2_2 - ((((epsivu2 * ewvsr - S * (epsilon))/ewusr2 - (0.5 +
    2 * (ewv/ewu2))) * depsi2 * ewvsr/ewusr2 + (0.5 * epsiu2 + 2 * (ewv/ewu2)) *
    pepsi2 - euv2 * sigx8_2) * sigx0_2 + 0.5 * (sigx7_2 * sigx2_2 - sigx8_2 *
    sigx0_2)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2) - sigx10_2 * (0.5 * ((2 *
    sigx4_2 + 2 * sigx4_1) * ewusr2) + 2 * (sigx10_2 * (1 - pwZ)))/((2 * sigx4_2 +
    2 * sigx4_1) * ewusr2)^2) * (1 - pwZ)), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - pwZ) * (2 * (((dmusig2 * ewvsr/(4 * (ewu2 * ewusr2)) - 2 * (ewu2 * pmusig2/(2 *
    ewu2)^2)) * ewv - ((0.5 * (epsiv2 * musig2) - 0.25) * dmusig2 * ewvsr/ewusr2 +
    (0.5 * epsiu2 + 2 * sigx6_2) * sigx9_2)) * sigx2_2 - ((2 * ((depsi2 * ewvsr/ewusr2 -
    pepsi2) * ewv/ewu2) - ((epsivu2 * sigx11_2 - 0.5) * depsi2 * ewvsr/ewusr2 +
    sigx12_2 * euv2)) * sigx0_2 + 0.5 * sigx13_2)) - 2 * (sigx10_2 * sigx16/(2 *
    sigx4_2 + 2 * sigx4_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx10_2 * (2 + 2 * ((1 - pwZ) * sigx17/(2 * sigx4_2 +
      2 * sigx4_1))) * dwZ/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2)), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (2 * ((((sigx9_1/2 + (pmusig1 - epsiv1 * dmusig1)/2) * ewv/ewu1 -
      (0.25 * (ewvsr/ewusr1) + 0.25 * (S * (epsilon)/ewvsr) - epsiv1^2 * musig1) *
        dmusig1) * sigx2_1 - ((2 * sigx12_1 + 2 * (pepsi1 - depsi1 * sigx11_1)) *
      ewv/ewu1 - (0.25 * (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr1) - epsivu1 *
      sigx11_1^2) * depsi1) * sigx0_1) * pwZ/ewusr1) + 2 * ((((sigx9_2/2 +
      (pmusig2 - epsiv2 * dmusig2)/2) * ewv/ewu2 - (0.25 * (ewvsr/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr) - epsiv2^2 * musig2) * dmusig2) * sigx2_2 -
      ((2 * sigx12_2 + 2 * (pepsi2 - depsi2 * sigx11_2)) * ewv/ewu2 - (0.25 *
        (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr2) - epsivu2 * sigx11_2^2) *
        depsi2) * sigx0_2) * (1 - pwZ)/ewusr2) - sigx16^2/(2 * sigx4_2 +
      2 * sigx4_1))/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx13_1/ewusr1) - (sigx16 * sigx17/(2 *
      sigx4_2 + 2 * sigx4_1) + 2 * (sigx13_2/ewusr2))) * dwZ/(2 * sigx4_2 +
      2 * sigx4_1), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx17 * dwZ/(2 * sigx4_2 + 2 * sigx4_1) +
      Wz) * sigx17 * dwZ/(2 * sigx4_2 + 2 * sigx4_1)), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cloglog specification class membership
chessmisfgenexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr <- exp(Wv/2)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr/ewusr1 + S * (epsilon)/ewvsr)
  musig2 <- (ewvsr/ewusr2 + S * (epsilon)/ewvsr)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewvsr/ewusr1) + S * (epsilon)/ewvsr)
  epsivu2 <- (2 * (ewvsr/ewusr2) + S * (epsilon)/ewvsr)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  dpuv1 <- (depsi1/ewvsr - 2 * (pepsi1/ewusr1))
  dpuv2 <- (depsi2/ewvsr - 2 * (pepsi2/ewusr2))
  euv1 <- (2 * (ewv/ewu1) + S * (epsilon)/ewusr1)
  euv2 <- (2 * (ewv/ewu2) + S * (epsilon)/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx0_1 <- exp(2 * (ewv/ewu1) + 2 * epsiu1)
  sigx0_2 <- exp(2 * (ewv/ewu2) + 2 * epsiu2)
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- ((sigx1_1 * sigx2_1 - dpuv1 * sigx0_1) * (1 - prZ)/ewusr1)
  sigx3_2 <- ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2) * prZ/ewusr2)
  sigx14_1 <- (sigx2_1 * pmusig1 - sigx0_1 * pepsi1)
  sigx14_2 <- (sigx2_2 * pmusig2 - sigx0_2 * pepsi2)
  sigx4_1 <- ((1 - prZ) * sigx14_1/ewusr1)
  sigx4_2 <- (prZ * sigx14_2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) * pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) * pmusig2)
  sigx8_1 <- (depsi1 * ewvsr/ewusr1 - euv1 * pepsi1)
  sigx8_2 <- (depsi2 * ewvsr/ewusr2 - euv2 * pepsi2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10_1 <- (sigx7_1 * sigx2_1 - (sigx8_1 * sigx0_1 + 0.5 * sigx14_1))
  sigx10_2 <- (sigx7_2 * sigx2_2 - (sigx8_2 * sigx0_2 + 0.5 * sigx14_2))
  sigx11_1 <- (ewvsr/ewusr1 - 0.5 * (S * (epsilon)/ewvsr))
  sigx11_2 <- (ewvsr/ewusr2 - 0.5 * (S * (epsilon)/ewvsr))
  sigx12_1 <- (2 * (ewv * pepsi1/ewu1) - depsi1 * sigx11_1)
  sigx12_2 <- (2 * (ewv * pepsi2/ewu2) - depsi2 * sigx11_2)
  sigx13_1 <- (sigx2_1 * sigx9_1 - sigx12_1 * sigx0_1)
  sigx13_2 <- (sigx2_2 * sigx9_2 - sigx12_2 * sigx0_2)
  sigx15 <- (pi * ((Wz)^2 + 1) * (2 * sigx4_2 + 2 * sigx4_1))
  sigx16 <- (2 * ((1 - prZ) * sigx13_1/ewusr1) + 2 * (prZ * sigx13_2/ewusr2))
  sigx17 <- (2 * (sigx14_1/ewusr1) - 2 * (sigx14_2/ewusr2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar, ncol = nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (2 * ((((musig1/ewvsr - 1/ewusr1) * dmusig1/ewvsr - sigx1_1/ewusr1) * sigx2_1 -
      ((epsivu1/ewvsr - 2/ewusr1) * depsi1/ewvsr - 2 * (dpuv1/ewusr1)) * sigx0_1) *
      (1 - prZ)/ewusr1) + 2 * ((((musig2/ewvsr - 1/ewusr2) * dmusig2/ewvsr -
      sigx1_2/ewusr2) * sigx2_2 - ((epsivu2/ewvsr - 2/ewusr2) * depsi2/ewvsr -
      2 * (dpuv2/ewusr2)) * sigx0_2) * prZ/ewusr2) - (2 * sigx3_1 + 2 * sigx3_2)^2/(2 *
      sigx4_2 + 2 * sigx4_1))/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * epsiu1 + 2 * sigx6_1) * pmusig1 -
      0.5 * sigx5_1)/ewusr1 + (0.5 * (musig1/ewusr1) - (0.5 * epsiu1 + 2 *
      sigx6_1)/ewvsr) * dmusig1) * sigx2_1 - (((epsivu1/ewusr1 - euv1/ewvsr) *
      depsi1 + (pepsi1 - 2 * sigx8_1)/ewusr1) * sigx0_1 + 0.5 * (sigx1_1 *
      sigx2_1 - dpuv1 * sigx0_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1) -
      sigx10_1 * (2 * sigx3_1 + 2 * sigx3_2) * ewusr1/((2 * sigx4_2 + 2 * sigx4_1) *
        ewusr1)^2) * (1 - prZ)), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
      pmusig2 - 0.5 * sigx5_2)/ewusr2 + (0.5 * (musig2/ewusr2) - (0.5 * epsiu2 +
      2 * sigx6_2)/ewvsr) * dmusig2) * sigx2_2 - (((epsivu2/ewusr2 - euv2/ewvsr) *
      depsi2 + (pepsi2 - 2 * sigx8_2)/ewusr2) * sigx0_2 + 0.5 * (sigx1_2 *
      sigx2_2 - dpuv2 * sigx0_2)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2) -
      sigx10_2 * (2 * sigx3_1 + 2 * sigx3_2) * ewusr2/((2 * sigx4_2 + 2 * sigx4_1) *
        ewusr2)^2) * prZ), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((dmusig1 * (ewv/(2 * ewu1) - (epsiv1 *
      musig1 + 0.5))/ewvsr - sigx9_1/ewusr1) * sigx2_1 - ((2 * (ewv/ewu1) -
      (epsivu1 * sigx11_1 + 0.5)) * depsi1/ewvsr - 2 * (sigx12_1/ewusr1)) *
      sigx0_1) * (1 - prZ)/ewusr1) + 2 * (((dmusig2 * (ewv/(2 * ewu2) - (epsiv2 *
      musig2 + 0.5))/ewvsr - sigx9_2/ewusr2) * sigx2_2 - ((2 * (ewv/ewu2) -
      (epsivu2 * sigx11_2 + 0.5)) * depsi2/ewvsr - 2 * (sigx12_2/ewusr2)) *
      sigx0_2) * prZ/ewusr2) - (2 * sigx3_1 + 2 * sigx3_2) * sigx16/(2 * sigx4_2 +
      2 * sigx4_1))/(2 * sigx4_2 + 2 * sigx4_1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * (2 * ((sigx1_1 *
    sigx2_1 - dpuv1 * sigx0_1)/ewusr1) - ((2 * sigx3_1 + 2 * sigx3_2) * sigx17/(2 *
    sigx4_1 + 2 * sigx4_2) + 2 * ((sigx1_2 * sigx2_2 - dpuv2 * sigx0_2)/ewusr2))) *
    prZ * ewz/(2 * sigx4_1 + 2 * sigx4_2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewvsr * musig1/ewusr1) -
      0.5) - 0.5 * (0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr/ewusr1 -
      (sigx7_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 - 8 * (ewu1^2/(2 *
        ewu1)^2)) * ewu1 * ewv/(2 * ewu1)^2) - 0.25 * epsiu1) * pmusig1)) *
      sigx2_1 - ((((epsivu1 * ewvsr - S * (epsilon))/ewusr1 - (0.5 + 2 * (ewv/ewu1))) *
      depsi1 * ewvsr/ewusr1 + (0.5 * epsiu1 + 2 * (ewv/ewu1)) * pepsi1 - euv1 *
      sigx8_1) * sigx0_1 + 0.5 * (sigx7_1 * sigx2_1 - sigx8_1 * sigx0_1)))/((2 *
      sigx4_2 + 2 * sigx4_1) * ewusr1) - sigx10_1 * (0.5 * ((2 * sigx4_2 +
      2 * sigx4_1) * ewusr1) + 2 * (sigx10_1 * (1 - prZ)))/((2 * sigx4_2 +
      2 * sigx4_1) * ewusr1)^2) * (1 - prZ)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (sigx10_1 * sigx10_2 * prZ * (1 - prZ) *
      ewusr2/(((2 * sigx4_2 + 2 * sigx4_1) * ewusr2)^2 * ewusr1))), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar + 2 *
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - prZ) * (2 * (((dmusig1 * ewvsr/(4 * (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 *
    ewu1)^2)) * ewv - ((0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 * ewvsr/ewusr1 +
    (0.5 * epsiu1 + 2 * sigx6_1) * sigx9_1)) * sigx2_1 - ((2 * ((depsi1 * ewvsr/ewusr1 -
    pepsi1) * ewv/ewu1) - ((epsivu1 * sigx11_1 - 0.5) * depsi1 * ewvsr/ewusr1 +
    sigx12_1 * euv1)) * sigx0_1 + 0.5 * sigx13_1)) - 2 * (sigx10_1 * sigx16/(2 *
    sigx4_2 + 2 * sigx4_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr1), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx10_1 * (2 - 2 * ((1 - prZ) * sigx17/(2 * sigx4_1 + 2 * sigx4_2))) * prZ *
    ewz/((2 * sigx4_1 + 2 * sigx4_2) * ewusr1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((((0.5 *
    (0.5 * (ewvsr * musig2/ewusr2) - 0.5) - 0.5 * (0.5 * epsiu2 + 2 * sigx6_2)) *
    dmusig2 * ewvsr/ewusr2 - (sigx7_2 * (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 -
    8 * (ewu2^2/(2 * ewu2)^2)) * ewu2 * ewv/(2 * ewu2)^2) - 0.25 * epsiu2) *
    pmusig2)) * sigx2_2 - ((((epsivu2 * ewvsr - S * (epsilon))/ewusr2 - (0.5 +
    2 * (ewv/ewu2))) * depsi2 * ewvsr/ewusr2 + (0.5 * epsiu2 + 2 * (ewv/ewu2)) *
    pepsi2 - euv2 * sigx8_2) * sigx0_2 + 0.5 * (sigx7_2 * sigx2_2 - sigx8_2 *
    sigx0_2)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2) - sigx10_2 * (0.5 * ((2 *
    sigx4_2 + 2 * sigx4_1) * ewusr2) + 2 * (sigx10_2 * prZ))/((2 * sigx4_2 +
    2 * sigx4_1) * ewusr2)^2) * prZ), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    prZ * (2 * (((dmusig2 * ewvsr/(4 * (ewu2 * ewusr2)) - 2 * (ewu2 * pmusig2/(2 *
    ewu2)^2)) * ewv - ((0.5 * (epsiv2 * musig2) - 0.25) * dmusig2 * ewvsr/ewusr2 +
    (0.5 * epsiu2 + 2 * sigx6_2) * sigx9_2)) * sigx2_2 - ((2 * ((depsi2 * ewvsr/ewusr2 -
    pepsi2) * ewv/ewu2) - ((epsivu2 * sigx11_2 - 0.5) * depsi2 * ewvsr/ewusr2 +
    sigx12_2 * euv2)) * sigx0_2 + 0.5 * sigx13_2)) - 2 * (sigx10_2 * sigx16/(2 *
    sigx4_2 + 2 * sigx4_1)))/((2 * sigx4_2 + 2 * sigx4_1) * ewusr2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx10_2 * (2 + 2 * (sigx17 * prZ/(2 * sigx4_1 + 2 * sigx4_2))) *
      prZ * ewz/((2 * sigx4_1 + 2 * sigx4_2) * ewusr2)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (2 * ((((sigx9_1/2 + (pmusig1 - epsiv1 * dmusig1)/2) * ewv/ewu1 -
      (0.25 * (ewvsr/ewusr1) + 0.25 * (S * (epsilon)/ewvsr) - epsiv1^2 * musig1) *
        dmusig1) * sigx2_1 - ((2 * sigx12_1 + 2 * (pepsi1 - depsi1 * sigx11_1)) *
      ewv/ewu1 - (0.25 * (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr1) - epsivu1 *
      sigx11_1^2) * depsi1) * sigx0_1) * (1 - prZ)/ewusr1) + 2 * ((((sigx9_2/2 +
      (pmusig2 - epsiv2 * dmusig2)/2) * ewv/ewu2 - (0.25 * (ewvsr/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr) - epsiv2^2 * musig2) * dmusig2) * sigx2_2 -
      ((2 * sigx12_2 + 2 * (pepsi2 - depsi2 * sigx11_2)) * ewv/ewu2 - (0.25 *
        (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr2) - epsivu2 * sigx11_2^2) *
        depsi2) * sigx0_2) * prZ/ewusr2) - sigx16^2/(2 * sigx4_2 + 2 * sigx4_1))/(2 *
      sigx4_2 + 2 * sigx4_1), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar), (nXvar + 2 *
    nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx13_1/ewusr1) - (sigx16 * sigx17/(2 *
      sigx4_1 + 2 * sigx4_2) + 2 * (sigx13_2/ewusr2))) * prZ * ewz/(2 * sigx4_1 +
      2 * sigx4_2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx17 * prZ/(2 * sigx4_1 + 2 * sigx4_2) +
      1) * ewz) * sigx17 * prZ * ewz/(2 * sigx4_1 + 2 * sigx4_2), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for misf generalized genexponential-normal distribution
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
misfgenexponormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfgenexponormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfgenexponormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfgenexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfgenexponormlike_logit,
    grad = cgradmisfgenexponormlike_logit, hess = chessmisfgenexponormlike_logit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfgenexponormlike_logit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfgenexponormlike_logit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfgenexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfgenexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfgenexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfgenexponormlike_logit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfgenexponormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfgenexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfgenexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfgenexponormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

# cauchit specification class membership
misfgenexponormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfgenexponormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfgenexponormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfgenexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfgenexponormlike_cauchit,
    grad = cgradmisfgenexponormlike_cauchit, hess = chessmisfgenexponormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfgenexponormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfgenexponormlike_cauchit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfgenexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfgenexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfgenexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfgenexponormlike_cauchit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfgenexponormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfgenexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfgenexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfgenexponormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

# probit specification class membership
misfgenexponormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfgenexponormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfgenexponormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfgenexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfgenexponormlike_probit,
    grad = cgradmisfgenexponormlike_probit, hess = chessmisfgenexponormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfgenexponormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfgenexponormlike_probit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfgenexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfgenexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfgenexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfgenexponormlike_probit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfgenexponormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfgenexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfgenexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfgenexponormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

# cloglog specification class membership
misfgenexponormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar,
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmisfgenexponormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmisfgenexponormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfgenexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = printInfo, maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmisfgenexponormlike_cloglog,
    grad = cgradmisfgenexponormlike_cloglog, hess = chessmisfgenexponormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmisfgenexponormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmisfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmisfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfgenexponormlike_cloglog(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmisfgenexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfgenexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chessmisfgenexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfgenexponormlike_cloglog(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessmisfgenexponormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfgenexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfgenexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfgenexponormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initGenExpo = initGenExpo))
}

# Conditional efficiencies estimation ----------
#' efficiencies for misf genexpo-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# logit specification class membership
cmisfgenexponormeff_logit <- function(object, level) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) - exp(B1) * (dnorm(b1) +
    b1 * pnorm(b1)))/(exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv/2) * (exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) - exp(B2) * (dnorm(b2) +
    b2 * pnorm(b2)))/(exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A1) * exp(-a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 - exp(Wv/2)) -
      exp(B1) * exp(-b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 - exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (exp(A2) * exp(-a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 - exp(Wv/2)) -
      exp(B2) * exp(-b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 - exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A1) * exp(a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 +
      exp(Wv/2)) - exp(B1) * exp(b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 + exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (exp(A2) * exp(a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 +
      exp(Wv/2)) - exp(B2) * exp(b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 + exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
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
cmisfgenexponormeff_cauchit <- function(object, level) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) - exp(B1) * (dnorm(b1) +
    b1 * pnorm(b1)))/(exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv/2) * (exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) - exp(B2) * (dnorm(b2) +
    b2 * pnorm(b2)))/(exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A1) * exp(-a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 - exp(Wv/2)) -
      exp(B1) * exp(-b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 - exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (exp(A2) * exp(-a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 - exp(Wv/2)) -
      exp(B2) * exp(-b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 - exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A1) * exp(a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 +
      exp(Wv/2)) - exp(B1) * exp(b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 + exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (exp(A2) * exp(a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 +
      exp(Wv/2)) - exp(B2) * exp(b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 + exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
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
cmisfgenexponormeff_probit <- function(object, level) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) - exp(B1) * (dnorm(b1) +
    b1 * pnorm(b1)))/(exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv/2) * (exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) - exp(B2) * (dnorm(b2) +
    b2 * pnorm(b2)))/(exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A1) * exp(-a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 - exp(Wv/2)) -
      exp(B1) * exp(-b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 - exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (exp(A2) * exp(-a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 - exp(Wv/2)) -
      exp(B2) * exp(-b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 - exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A1) * exp(a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 +
      exp(Wv/2)) - exp(B1) * exp(b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 + exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (exp(A2) * exp(a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 +
      exp(Wv/2)) - exp(B2) * exp(b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 + exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
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
cmisfgenexponormeff_cloglog <- function(object, level) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) - exp(B1) * (dnorm(b1) +
    b1 * pnorm(b1)))/(exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv/2) * (exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) - exp(B2) * (dnorm(b2) +
    b2 * pnorm(b2)))/(exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A1) * exp(-a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 - exp(Wv/2)) -
      exp(B1) * exp(-b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 - exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (exp(A2) * exp(-a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 - exp(Wv/2)) -
      exp(B2) * exp(-b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 - exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A1) * exp(a1 * exp(Wv/2) + exp(Wv)/2) * pnorm(a1 +
      exp(Wv/2)) - exp(B1) * exp(b1 * exp(Wv/2) + exp(Wv)/2) * pnorm(b1 + exp(Wv/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (exp(A2) * exp(a2 * exp(Wv/2) + exp(Wv)/2) * pnorm(a2 +
      exp(Wv/2)) - exp(B2) * exp(b2 * exp(Wv/2) + exp(Wv)/2) * pnorm(b2 + exp(Wv/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
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
#' marginal impact on efficiencies for misf genexpo-normal distribution
#' @param object object of class sfacross
#' @noRd
# logit specification class membership
cmisfmarggenexponorm_Eu_logit <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu1/2),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu2/2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarggenexponorm_Vu_logit <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cauchit specification class membership
cmisfmarggenexponorm_Eu_cauchit <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu1/2),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu2/2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarggenexponorm_Vu_cauchit <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# probit specification class membership
cmisfmarggenexponorm_Eu_probit <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu1/2),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu2/2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarggenexponorm_Vu_probit <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# cloglog specification class membership
cmisfmarggenexponorm_Eu_cloglog <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu1/2),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 3/4, nrow = 1), matrix(exp(Wu2/2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmisfmarggenexponorm_Vu_cloglog <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 5/4, nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
