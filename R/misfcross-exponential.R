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
# Convolution: exponential - normal                                            #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for misf exponential-normal distribution
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
cmisfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# cauchit specification class membership
cmisfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# probit specification class membership
cmisfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# cloglog specification class membership
cmisfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for misf exponential-normal distribution
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
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
cstmisfexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + exponential - normal distributions...\n")
  initExpo <- maxLik(logLik = cexponormlike, start = cstexponorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradexponormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = printInfo, reltol = tol), nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]),
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar)
  Esti <- initExpo$estimate
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("MI_", colnames(Zvar)))
  names(initExpo$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initExpo = initExpo))
}

# Gradient of the likelihood function ----------
#' gradient for misf exponential-normal distribution
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
cgradmisfexponormlike_logit <- function(parm, nXvar, nuZUvar,
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
  ewv_h <- exp(Wv/2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv <- exp(Wv)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv_h/ewu1_h + S * (epsilon)/ewv_h)
  musig2 <- (ewv_h/ewu2_h + S * (epsilon)/ewv_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx1_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx2_1 <- (dmusig1/ewv_h - pmusig1/ewu1_h)
  sigx2_2 <- (dmusig2/ewv_h - pmusig2/ewu2_h)
  wzu1 <- (wzdeno * ewu1_h)
  wzu2 <- (wzdeno * ewu2_h)
  sigx3 <- (prC * sigx2_2 * sigx1_2/ewu2_h + sigx2_1 * sigx1_1 *
    ewz/wzu1)
  sigx4 <- (prC * sigx1_2 * pmusig2/ewu2_h + sigx1_1 * ewz *
    pmusig1/wzu1)
  duv1 <- (dmusig1 * ewv_h/ewu1_h)
  duv2 <- (dmusig2 * ewv_h/ewu2_h)
  ueps1 <- (S * (epsilon)/ewu1_h)
  ueps2 <- (S * (epsilon)/ewu2_h)
  usq1 <- (ewu1 * ewv/(2 * ewu1)^2)
  usq2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx5_1 <- (0.5 * duv1 - (0.5 * ueps1 + 2 * usq1) * pmusig1)
  sigx5_2 <- (0.5 * duv2 - (0.5 + 0.5 * ueps2 + 2 * usq2) *
    pmusig2)
  sigx6_1 <- (wzdeno * ewu1_h * pmusig1/wzu1^2)
  sigx7_1 <- (sigx5_1/wzu1 - 0.5 * sigx6_1)
  vu1 <- (ewv_h/ewu1_h)
  vu2 <- (ewv_h/ewu2_h)
  epsiv1 <- (S * (epsilon)/ewv_h)
  epsiv2 <- (S * (epsilon)/ewv_h)
  sigx8_1 <- (ewv * pmusig1/(2 * ewu1) - (0.5 * vu1 - 0.5 *
    epsiv1) * dmusig1)
  sigx8_2 <- (ewv * pmusig2/(2 * ewu2) - (0.5 * vu2 - 0.5 *
    epsiv2) * dmusig2)
  wzsq <- (1/wzu1 - ewu1_h * ewz/wzu1^2)
  sigx9 <- (wzsq * sigx1_1 * pmusig1 - prC * sigx1_2 * pmusig2/wzu2)
  s4u1 <- (sigx4 * wzdeno * ewu1_h)
  s4u2 <- (sigx4 * ewu2_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx7_1 *
    sigx1_1 * ewz/sigx4, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = sigx5_2 * prC * sigx1_2/s4u2, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = sigx1_1 * sigx8_1 * ewz/s4u1 + prC *
      sigx1_2 * sigx8_2/s4u2, FUN = "*"), sweep(Zvar, MARGIN = 1,
    STATS = sigx9 * ewz/sigx4, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cauchit specification class membership
cgradmisfexponormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- (ewz2 * sigx1_2 * sigx2_2/ewusr2 + ewz1 * sigx1_1 *
    sigx2_1/ewusr1)
  sigx4 <- (ewz2 * sigx2_2 * pmusig2/ewusr2 + ewz1 * sigx2_1 *
    pmusig1/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  sigx11 <- (ewz2 * sigx2_2 * sigx9_2/ewusr2 + ewz1 * sigx2_1 *
    sigx9_1/ewusr1)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx8_1 *
    ewz1 * sigx2_1/(sigx4 * ewusr1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = ewz2 * sigx8_2 * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx11/sigx4,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx10/(pi *
    sigx4 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# probit specification class membership
cgradmisfexponormlike_probit <- function(parm, nXvar, nuZUvar,
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
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - pwZ) * sigx1_2 * sigx2_2/ewusr2 + sigx1_1 *
    sigx2_1 * pwZ/ewusr1)
  sigx4 <- ((1 - pwZ) * sigx2_2 * pmusig2/ewusr2 + sigx2_1 *
    pmusig1 * pwZ/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  sigx11 <- ((1 - pwZ) * sigx2_2 * sigx9_2/ewusr2 + sigx2_1 *
    sigx9_1 * pwZ/ewusr1)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx8_1 *
    pwZ * sigx2_1/(sigx4 * ewusr1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = (1 - pwZ) * sigx8_2 * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx11/sigx4,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = dwZ * sigx10/sigx4,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# cloglog specification class membership
cgradmisfexponormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - prZ) * sigx1_1 * sigx2_1/ewusr1 + sigx1_2 *
    prZ * sigx2_2/ewusr2)
  sigx4 <- ((1 - prZ) * sigx2_1 * pmusig1/ewusr1 + prZ * sigx2_2 *
    pmusig2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  sigx11 <- ((1 - prZ) * sigx2_1 * sigx9_1/ewusr1 + prZ * sigx2_2 *
    sigx9_2/ewusr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx8_1 *
    (1 - prZ) * sigx2_1/(sigx4 * ewusr1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = prZ * sigx8_2 * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx11/sigx4,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = prZ * sigx10 *
    ewz/sigx4, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for misf exponential-normal distribution
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
chessmisfexponormlike_logit <- function(parm, nXvar, nuZUvar,
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
  ewv_h <- exp(Wv/2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv <- exp(Wv)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv_h/ewu1_h + S * (epsilon)/ewv_h)
  musig2 <- (ewv_h/ewu2_h + S * (epsilon)/ewv_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx1_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx2_1 <- (dmusig1/ewv_h - pmusig1/ewu1_h)
  sigx2_2 <- (dmusig2/ewv_h - pmusig2/ewu2_h)
  wzu1 <- (wzdeno * ewu1_h)
  wzu2 <- (wzdeno * ewu2_h)
  sigx3 <- (prC * sigx2_2 * sigx1_2/ewu2_h + sigx2_1 * sigx1_1 *
    ewz/wzu1)
  sigx4 <- (prC * sigx1_2 * pmusig2/ewu2_h + sigx1_1 * ewz *
    pmusig1/wzu1)
  duv1 <- (dmusig1 * ewv_h/ewu1_h)
  duv2 <- (dmusig2 * ewv_h/ewu2_h)
  ueps1 <- (S * (epsilon)/ewu1_h)
  ueps2 <- (S * (epsilon)/ewu2_h)
  usq1 <- (ewu1 * ewv/(2 * ewu1)^2)
  usq2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx5_1 <- (0.5 * duv1 - (0.5 * ueps1 + 2 * usq1) * pmusig1)
  sigx5_2 <- (0.5 * duv2 - (0.5 + 0.5 * ueps2 + 2 * usq2) *
    pmusig2)
  sigx6_1 <- (wzdeno * ewu1_h * pmusig1/wzu1^2)
  sigx7_1 <- (sigx5_1/wzu1 - 0.5 * sigx6_1)
  vu1 <- (ewv_h/ewu1_h)
  vu2 <- (ewv_h/ewu2_h)
  epsiv1 <- (S * (epsilon)/ewv_h)
  epsiv2 <- (S * (epsilon)/ewv_h)
  sigx8_1 <- (ewv * pmusig1/(2 * ewu1) - (0.5 * vu1 - 0.5 *
    epsiv1) * dmusig1)
  sigx8_2 <- (ewv * pmusig2/(2 * ewu2) - (0.5 * vu2 - 0.5 *
    epsiv2) * dmusig2)
  wzsq <- (1/wzu1 - ewu1_h * ewz/wzu1^2)
  sigx9 <- (wzsq * sigx1_1 * pmusig1 - prC * sigx1_2 * pmusig2/wzu2)
  s4u2 <- (sigx4 * ewu2_h)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (((musig1/ewv_h - 1/ewu1_h) * dmusig1/ewv_h -
      sigx2_1/ewu1_h) * sigx1_1 * ewz/wzu1 + ((musig2/ewv_h -
      1/ewu2_h) * dmusig2/ewv_h - sigx2_2/ewu2_h) * prC *
      sigx1_2/ewu2_h - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((0.5 + 0.5 * ueps1 +
      2 * usq1) * pmusig1 - 0.5 * duv1)/ewu1_h + (0.5 *
      (musig1/ewu1_h) - (0.5 * ueps1 + 2 * usq1)/ewv_h) *
      dmusig1)/wzdeno + 0.5 * sigx6_1)/ewu1_h - (sigx7_1 *
      sigx3/sigx4 + 0.5 * (wzdeno * dmusig1 * ewu1_h/(wzu1^2 *
      ewv_h)))) * sigx1_1 * ewz/sigx4, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ueps2 + 1 +
      2 * usq2) * pmusig2 - 0.5 * duv2)/ewu2_h + (0.5 *
      (musig2/ewu2_h) - (0.5 + 0.5 * ueps2 + 2 * usq2)/ewv_h) *
      dmusig2)/s4u2 - sigx3 * sigx5_2 * ewu2_h/s4u2^2) *
      prC * sigx1_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (prC * (dmusig2 * (ewv/(2 * ewu2) - ((0.5 * vu2 -
    0.5 * epsiv1) * musig2 + 0.5))/ewv_h - sigx8_2/ewu2_h) *
    sigx1_2/ewu2_h + (dmusig1 * (ewv/(2 * ewu1) - ((0.5 *
    vu1 - 0.5 * epsiv1) * musig1 + 0.5))/ewv_h - sigx8_1/ewu1_h) *
    sigx1_1 * ewz/wzu1 - sigx3 * (prC * sigx1_2 * sigx8_2/ewu2_h +
    sigx1_1 * sigx8_1 * ewz/wzu1)/sigx4)/sigx4, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (wzsq * sigx2_1 * sigx1_1 -
      (sigx3 * sigx9/sigx4 + prC * sigx2_2 * sigx1_2/wzu2)) *
      ewz/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewv_h * musig1/ewu1_h) - 0.5) - 0.5 *
      (0.5 * ueps1 + 2 * usq1)) * dmusig1 * ewv_h/ewu1_h -
      (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv/(2 *
        ewu1)^2) - 0.25 * ueps1) * pmusig1)/wzu1 - ((sigx7_1 *
      sigx1_1 * ewz/sigx4 + 0.5 * ueps1 + 2 * usq1) * sigx7_1 +
      (0.5 * ((0.5 - wzdeno^2 * ewu1_h^2/wzu1^2) * ewu1_h *
        pmusig1 + 0.5 * (dmusig1 * ewv_h)) + 0.5 * (sigx5_1 *
        ewu1_h)) * wzdeno/wzu1^2)) * sigx1_1 * ewz/sigx4,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx7_1 * sigx5_2 * prC * ewu2_h *
      sigx1_1 * sigx1_2 * ewz/s4u2^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((dmusig1 * ewv_h/(4 * (ewu1 *
      ewu1_h)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) * ewv -
      ((0.5 * ((0.5 * vu1 - 0.5 * epsiv1) * musig1) - 0.25) *
        dmusig1 * ewv_h/ewu1_h + (0.5 * ueps1 + 2 * usq1) *
        sigx8_1))/wzu1 - (sigx7_1 * (prC * sigx1_2 *
      sigx8_2/ewu2_h + sigx1_1 * sigx8_1 * ewz/wzu1)/sigx4 +
      0.5 * (wzdeno * ewu1_h * sigx8_1/wzu1^2))) * sigx1_1 *
      ewz/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (0.5 * (wzsq * dmusig1 *
      ewv_h/ewu1_h) - ((((0.5 - wzdeno^2 * ewu1_h^2/wzu1^2) *
      ewz + 0.5 * wzdeno) * ewu1_h/wzu1^2 + (0.5 * ueps1 +
      2 * usq1) * wzsq) * pmusig1 + sigx7_1 * sigx9 * ewz/sigx4)) *
      sigx1_1 * ewz/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewv_h *
      musig2/ewu2_h) - 0.5) - 0.5 * (0.5 + 0.5 * ueps2 +
      2 * usq2)) * dmusig2 * ewv_h/ewu2_h - (sigx5_2 *
      (0.5 * ueps2 + 2 * usq2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv/(2 * ewu2)^2) - 0.25 * ueps2) *
      pmusig2))/s4u2 - (sigx5_2 * prC * sigx1_2 + 0.5 *
      s4u2) * sigx5_2/s4u2^2) * prC * sigx1_2, FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewv_h/(4 * (ewu2 *
      ewu2_h)) - 2 * (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv -
      ((prC * sigx1_2 * sigx8_2/ewu2_h + sigx1_1 * sigx8_1 *
        ewz/wzu1) * sigx5_2/sigx4 + (0.5 * ((0.5 * vu2 -
        0.5 * epsiv1) * musig2) - 0.25) * dmusig2 * ewv_h/ewu2_h +
        (0.5 + 0.5 * ueps2 + 2 * usq2) * sigx8_2)) *
      prC * sigx1_2/s4u2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((sigx9 * sigx5_2/sigx4 + (0.5 * duv2 - (0.5 * ueps2 +
      2 * usq2) * pmusig2)/wzdeno)/ewu2_h - 0.5 * (wzdeno *
      ewu2_h * pmusig2/wzu2^2)) * prC * sigx1_2 * ewz/sigx4),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx8_1/2 + (pmusig1 -
      (0.5 * vu1 - 0.5 * epsiv1) * dmusig1)/2) * ewv/ewu1 -
      (0.25 * vu1 + 0.25 * epsiv1 - (0.5 * vu1 - 0.5 *
        epsiv1)^2 * musig1) * dmusig1) * sigx1_1 * ewz/wzu1 +
      ((sigx8_2/2 + (pmusig2 - (0.5 * vu2 - 0.5 * epsiv1) *
        dmusig2)/2) * ewv/ewu2 - (0.25 * vu2 + 0.25 *
        epsiv1 - (0.5 * vu2 - 0.5 * epsiv1)^2 * musig2) *
        dmusig2) * prC * sigx1_2/ewu2_h - (prC * sigx1_2 *
      sigx8_2/ewu2_h + sigx1_1 * sigx8_1 * ewz/wzu1)^2/sigx4)/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (wzsq * sigx1_1 * sigx8_1 - ((prC * sigx1_2 *
      sigx8_2/ewu2_h + sigx1_1 * sigx8_1 * ewz/wzu1) *
      sigx9/sigx4 + prC * sigx1_2 * sigx8_2/wzu2)) * ewz/sigx4,
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewu2_h) +
      ewu2_h/wzu2^2) * sigx1_2 * pmusig2 - (sigx9^2/sigx4 +
      (2 - 2 * (wzdeno * ewu1_h^2 * ewz/wzu1^2)) * ewu1_h *
        sigx1_1 * pmusig1/wzu1^2)) * ewz + wzsq * sigx1_1 *
      pmusig1 - prC * sigx1_2 * pmusig2/wzu2) * ewz/sigx4,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cauchit specification class membership
chessmisfexponormlike_cauchit <- function(parm, nXvar, nuZUvar,
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
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- (ewz2 * sigx1_2 * sigx2_2/ewusr2 + ewz1 * sigx1_1 *
    sigx2_1/ewusr1)
  sigx4 <- (ewz2 * sigx2_2 * pmusig2/ewusr2 + ewz1 * sigx2_1 *
    pmusig1/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  sigx11 <- (ewz2 * sigx2_2 * sigx9_2/ewusr2 + ewz1 * sigx2_1 *
    sigx9_1/ewusr1)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr - 1/ewusr1) * dmusig1/ewvsr -
      sigx1_1/ewusr1) * ewz1 * sigx2_1/ewusr1 + ((musig2/ewvsr -
      1/ewusr2) * dmusig2/ewvsr - sigx1_2/ewusr2) * ewz2 *
      sigx2_2/ewusr2 - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu1 + 1 +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 + 0.5 * epsiu1 + 2 *
        sigx6_1)/ewvsr) * dmusig1)/(sigx4 * ewusr1) -
      sigx3 * sigx8_1 * ewusr1/(sigx4 * ewusr1)^2) * ewz1 *
      sigx2_1, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu2 + 1 +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 + 0.5 * epsiu2 + 2 *
        sigx6_2)/ewvsr) * dmusig2)/(sigx4 * ewusr2) -
      sigx3 * sigx8_2 * ewusr2/(sigx4 * ewusr2)^2) * ewz2 *
      sigx2_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (ewz2 * (dmusig2 * (ewv/(2 * ewu2) - (epsiv2 * musig2 +
    0.5))/ewvsr - sigx9_2/ewusr2) * sigx2_2/ewusr2 + ewz1 *
    (dmusig1 * (ewv/(2 * ewu1) - (epsiv1 * musig1 + 0.5))/ewvsr -
      sigx9_1/ewusr1) * sigx2_1/ewusr1 - sigx3 * sigx11/sigx4)/sigx4,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx1_1 * sigx2_1/ewusr1 -
      sigx1_2 * sigx2_2/ewusr2)/(pi * sigx4 * ((Wz)^2 +
      1)) - pi * sigx3 * ((Wz)^2 + 1) * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr * musig1/ewusr1) - 0.5) - 0.5 *
      (0.5 + 0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr/ewusr1 -
      (sigx8_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv/(2 *
        ewu1)^2) - 0.25 * epsiu1) * pmusig1))/(sigx4 *
      ewusr1) - (sigx8_1 * ewz1 * sigx2_1 + 0.5 * (sigx4 *
      ewusr1)) * sigx8_1/(sigx4 * ewusr1)^2) * ewz1 * sigx2_1,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (ewz2 * sigx8_1 * sigx8_2 * ewz1 * ewusr2 *
      sigx2_1 * sigx2_2/((sigx4 * ewusr2)^2 * ewusr1)),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig1 * ewvsr/(4 * (ewu1 *
      ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) * ewv -
      (sigx11 * sigx8_1/sigx4 + (0.5 * (epsiv1 * musig1) -
        0.25) * dmusig1 * ewvsr/ewusr1 + (0.5 + 0.5 *
        epsiu1 + 2 * sigx6_1) * sigx9_1)) * ewz1 * sigx2_1/(sigx4 *
      ewusr1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx8_1 * (1/(pi * sigx4 *
      ((Wz)^2 + 1)) - pi * ((Wz)^2 + 1) * ewz1 * sigx10/(pi *
      sigx4 * ((Wz)^2 + 1))^2) * sigx2_1/ewusr1, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 + 0.5 * epsiu2 +
      2 * sigx6_2)) * dmusig2 * ewvsr/ewusr2 - (sigx8_2 *
      (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv/(2 * ewu2)^2) - 0.25 * epsiu2) *
      pmusig2))/(sigx4 * ewusr2) - (ewz2 * sigx8_2 * sigx2_2 +
      0.5 * (sigx4 * ewusr2)) * sigx8_2/(sigx4 * ewusr2)^2) *
      ewz2 * sigx2_2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewvsr/(4 * (ewu2 *
      ewusr2)) - 2 * (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv -
      (sigx11 * sigx8_2/sigx4 + (0.5 * (epsiv2 * musig2) -
        0.25) * dmusig2 * ewvsr/ewusr2 + (0.5 + 0.5 *
        epsiu2 + 2 * sigx6_2) * sigx9_2)) * ewz2 * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx8_2 * (1/(pi * sigx4 * ((Wz)^2 + 1)) + pi * ((Wz)^2 +
      1) * ewz2 * sigx10/(pi * sigx4 * ((Wz)^2 + 1))^2) *
      sigx2_2/ewusr2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv/ewu1 - (0.25 * (ewvsr/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr) - epsiv1^2 * musig1) *
      dmusig1) * ewz1 * sigx2_1/ewusr1 + ((sigx9_2/2 +
      (pmusig2 - epsiv2 * dmusig2)/2) * ewv/ewu2 - (0.25 *
      (ewvsr/ewusr2) + 0.25 * (S * (epsilon)/ewvsr) - epsiv2^2 *
      musig2) * dmusig2) * ewz2 * sigx2_2/ewusr2 - sigx11^2/sigx4)/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx2_1 * sigx9_1/ewusr1 - sigx2_2 *
      sigx9_2/ewusr2)/(pi * sigx4 * ((Wz)^2 + 1)) - pi *
      sigx11 * ((Wz)^2 + 1) * sigx10/(pi * sigx4 * ((Wz)^2 +
      1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx4) +
      sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2) *
      sigx10/(pi * sigx4 * ((Wz)^2 + 1))^2), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# probit specification class membership
chessmisfexponormlike_probit <- function(parm, nXvar, nuZUvar,
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
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - pwZ) * sigx1_2 * sigx2_2/ewusr2 + sigx1_1 *
    sigx2_1 * pwZ/ewusr1)
  sigx4 <- ((1 - pwZ) * sigx2_2 * pmusig2/ewusr2 + sigx2_1 *
    pmusig1 * pwZ/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  sigx11 <- ((1 - pwZ) * sigx2_2 * sigx9_2/ewusr2 + sigx2_1 *
    sigx9_1 * pwZ/ewusr1)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr - 1/ewusr1) * dmusig1/ewvsr -
      sigx1_1/ewusr1) * pwZ * sigx2_1/ewusr1 + ((musig2/ewvsr -
      1/ewusr2) * dmusig2/ewvsr - sigx1_2/ewusr2) * (1 -
      pwZ) * sigx2_2/ewusr2 - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu1 + 1 +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 + 0.5 * epsiu1 + 2 *
        sigx6_1)/ewvsr) * dmusig1)/(sigx4 * ewusr1) -
      sigx3 * sigx8_1 * ewusr1/(sigx4 * ewusr1)^2) * pwZ *
      sigx2_1, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu2 + 1 +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 + 0.5 * epsiu2 + 2 *
        sigx6_2)/ewvsr) * dmusig2)/(sigx4 * ewusr2) -
      sigx3 * sigx8_2 * ewusr2/(sigx4 * ewusr2)^2) * (1 -
      pwZ) * sigx2_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((1 - pwZ) * (dmusig2 * (ewv/(2 * ewu2) - (epsiv2 *
    musig2 + 0.5))/ewvsr - sigx9_2/ewusr2) * sigx2_2/ewusr2 +
    pwZ * (dmusig1 * (ewv/(2 * ewu1) - (epsiv1 * musig1 +
      0.5))/ewvsr - sigx9_1/ewusr1) * sigx2_1/ewusr1 -
    sigx3 * sigx11/sigx4)/sigx4, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1 * sigx2_1/ewusr1 -
      (sigx3 * sigx10/sigx4 + sigx1_2 * sigx2_2/ewusr2)) *
      dwZ/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr * musig1/ewusr1) - 0.5) - 0.5 *
      (0.5 + 0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr/ewusr1 -
      (sigx8_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv/(2 *
        ewu1)^2) - 0.25 * epsiu1) * pmusig1))/(sigx4 *
      ewusr1) - (sigx8_1 * pwZ * sigx2_1 + 0.5 * (sigx4 *
      ewusr1)) * sigx8_1/(sigx4 * ewusr1)^2) * pwZ * sigx2_1,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((1 - pwZ) * sigx8_1 * sigx8_2 * pwZ *
      ewusr2 * sigx2_1 * sigx2_2/((sigx4 * ewusr2)^2 *
      ewusr1)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig1 * ewvsr/(4 * (ewu1 *
      ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) * ewv -
      (sigx11 * sigx8_1/sigx4 + (0.5 * (epsiv1 * musig1) -
        0.25) * dmusig1 * ewvsr/ewusr1 + (0.5 + 0.5 *
        epsiu1 + 2 * sigx6_1) * sigx9_1)) * pwZ * sigx2_1/(sigx4 *
      ewusr1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx8_1 * (1 - sigx10 * pwZ/sigx4) *
      dwZ * sigx2_1/(sigx4 * ewusr1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 + 0.5 * epsiu2 +
      2 * sigx6_2)) * dmusig2 * ewvsr/ewusr2 - (sigx8_2 *
      (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv/(2 * ewu2)^2) - 0.25 * epsiu2) *
      pmusig2))/(sigx4 * ewusr2) - ((1 - pwZ) * sigx8_2 *
      sigx2_2 + 0.5 * (sigx4 * ewusr2)) * sigx8_2/(sigx4 *
      ewusr2)^2) * (1 - pwZ) * sigx2_2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewvsr/(4 * (ewu2 *
      ewusr2)) - 2 * (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv -
      (sigx11 * sigx8_2/sigx4 + (0.5 * (epsiv2 * musig2) -
        0.25) * dmusig2 * ewvsr/ewusr2 + (0.5 + 0.5 *
        epsiu2 + 2 * sigx6_2) * sigx9_2)) * (1 - pwZ) *
      sigx2_2/(sigx4 * ewusr2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((1 - pwZ) * sigx10/sigx4 + 1) * sigx8_2 * dwZ * sigx2_2/(sigx4 *
      ewusr2)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv/ewu1 - (0.25 * (ewvsr/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr) - epsiv1^2 * musig1) *
      dmusig1) * pwZ * sigx2_1/ewusr1 + ((sigx9_2/2 + (pmusig2 -
      epsiv2 * dmusig2)/2) * ewv/ewu2 - (0.25 * (ewvsr/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr) - epsiv2^2 * musig2) *
      dmusig2) * (1 - pwZ) * sigx2_2/ewusr2 - sigx11^2/sigx4)/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * dwZ * (sigx2_1 * sigx9_1/ewusr1 - (sigx11 *
      sigx10/sigx4 + sigx2_2 * sigx9_2/ewusr2))/sigx4,
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (dwZ * (dwZ * sigx10/sigx4 +
      Wz) * sigx10/sigx4), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# cloglog specification class membership
chessmisfexponormlike_cloglog <- function(parm, nXvar, nuZUvar,
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
  sigx1_1 <- (dmusig1/ewvsr - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - prZ) * sigx1_1 * sigx2_1/ewusr1 + sigx1_2 *
    prZ * sigx2_2/ewusr2)
  sigx4 <- ((1 - prZ) * sigx2_1 * pmusig1/ewusr1 + prZ * sigx2_2 *
    pmusig2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr/ewusr1) - 0.5 * (S * (epsilon)/ewvsr))
  epsiv2 <- (0.5 * (ewvsr/ewusr2) - 0.5 * (S * (epsilon)/ewvsr))
  sigx9_1 <- (ewv * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  sigx11 <- ((1 - prZ) * sigx2_1 * sigx9_1/ewusr1 + prZ * sigx2_2 *
    sigx9_2/ewusr2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + 2 * nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr - 1/ewusr1) * dmusig1/ewvsr -
      sigx1_1/ewusr1) * (1 - prZ) * sigx2_1/ewusr1 + ((musig2/ewvsr -
      1/ewusr2) * dmusig2/ewvsr - sigx1_2/ewusr2) * prZ *
      sigx2_2/ewusr2 - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu1 + 1 +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 + 0.5 * epsiu1 + 2 *
        sigx6_1)/ewvsr) * dmusig1)/(sigx4 * ewusr1) -
      sigx3 * sigx8_1 * ewusr1/(sigx4 * ewusr1)^2) * (1 -
      prZ) * sigx2_1, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu2 + 1 +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 + 0.5 * epsiu2 + 2 *
        sigx6_2)/ewvsr) * dmusig2)/(sigx4 * ewusr2) -
      sigx3 * sigx8_2 * ewusr2/(sigx4 * ewusr2)^2) * prZ *
      sigx2_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (prZ * (dmusig2 * (ewv/(2 * ewu2) - (epsiv2 * musig2 +
    0.5))/ewvsr - sigx9_2/ewusr2) * sigx2_2/ewusr2 + (1 -
    prZ) * (dmusig1 * (ewv/(2 * ewu1) - (epsiv1 * musig1 +
    0.5))/ewvsr - sigx9_1/ewusr1) * sigx2_1/ewusr1 - sigx3 *
    sigx11/sigx4)/sigx4, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1 * sigx2_1/ewusr1 -
      (sigx3 * sigx10/sigx4 + sigx1_2 * sigx2_2/ewusr2)) *
      prZ * ewz/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr * musig1/ewusr1) - 0.5) - 0.5 *
      (0.5 + 0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr/ewusr1 -
      (sigx8_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv/(2 *
        ewu1)^2) - 0.25 * epsiu1) * pmusig1))/(sigx4 *
      ewusr1) - (sigx8_1 * (1 - prZ) * sigx2_1 + 0.5 *
      (sigx4 * ewusr1)) * sigx8_1/(sigx4 * ewusr1)^2) *
    (1 - prZ) * sigx2_1, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (prZ * sigx8_1 * sigx8_2 * (1 - prZ) *
      ewusr2 * sigx2_1 * sigx2_2/((sigx4 * ewusr2)^2 *
      ewusr1)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig1 * ewvsr/(4 * (ewu1 *
      ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) * ewv -
      (sigx11 * sigx8_1/sigx4 + (0.5 * (epsiv1 * musig1) -
        0.25) * dmusig1 * ewvsr/ewusr1 + (0.5 + 0.5 *
        epsiu1 + 2 * sigx6_1) * sigx9_1)) * (1 - prZ) *
      sigx2_1/(sigx4 * ewusr1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx8_1 * (1 - (1 - prZ) *
      sigx10/sigx4) * prZ * sigx2_1 * ewz/(sigx4 * ewusr1),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 + 0.5 * epsiu2 +
      2 * sigx6_2)) * dmusig2 * ewvsr/ewusr2 - (sigx8_2 *
      (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv/(2 * ewu2)^2) - 0.25 * epsiu2) *
      pmusig2))/(sigx4 * ewusr2) - (prZ * sigx8_2 * sigx2_2 +
      0.5 * (sigx4 * ewusr2)) * sigx8_2/(sigx4 * ewusr2)^2) *
      prZ * sigx2_2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewvsr/(4 * (ewu2 *
      ewusr2)) - 2 * (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv -
      (sigx11 * sigx8_2/sigx4 + (0.5 * (epsiv2 * musig2) -
        0.25) * dmusig2 * ewvsr/ewusr2 + (0.5 + 0.5 *
        epsiu2 + 2 * sigx6_2) * sigx9_2)) * prZ * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx8_2 * (1 + prZ * sigx10/sigx4) * prZ * sigx2_2 *
      ewz/(sigx4 * ewusr2)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv/ewu1 - (0.25 * (ewvsr/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr) - epsiv1^2 * musig1) *
      dmusig1) * (1 - prZ) * sigx2_1/ewusr1 + ((sigx9_2/2 +
      (pmusig2 - epsiv2 * dmusig2)/2) * ewv/ewu2 - (0.25 *
      (ewvsr/ewusr2) + 0.25 * (S * (epsilon)/ewvsr) - epsiv2^2 *
      musig2) * dmusig2) * prZ * sigx2_2/ewusr2 - sigx11^2/sigx4)/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prZ * (sigx2_1 * sigx9_1/ewusr1 - (sigx11 *
      sigx10/sigx4 + sigx2_2 * sigx9_2/ewusr2)) * ewz/sigx4,
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar + nZHvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 + prZ * sigx10/sigx4) *
      ewz) * prZ * sigx10 * ewz/sigx4, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for misf exponential-normal distribution
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
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
# logit specification class membership
misfexponormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfexponormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfexponormlike_logit,
      grad = cgradmisfexponormlike_logit, hess = chessmisfexponormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfexponormlike_logit(mleObj$par,
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
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmisfexponormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfexponormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfexponormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

# cauchit specification class membership
misfexponormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfexponormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfexponormlike_cauchit,
      grad = cgradmisfexponormlike_cauchit, hess = chessmisfexponormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfexponormlike_cauchit(mleObj$par,
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
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmisfexponormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfexponormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfexponormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

# probit specification class membership
misfexponormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfexponormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfexponormlike_probit,
      grad = cgradmisfexponormlike_probit, hess = chessmisfexponormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfexponormlike_probit(mleObj$par,
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
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmisfexponormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfexponormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfexponormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

# cloglog specification class membership
misfexponormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmisfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmisfexponormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmisfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmisfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmisfexponormlike_cloglog,
      grad = cgradmisfexponormlike_cloglog, hess = chessmisfexponormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmisfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmisfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmisfexponormlike_cloglog(mleObj$par,
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
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmisfexponormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmisfexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmisfexponormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmisfexponormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

# Conditional efficiencies estimation ----------
#' efficiencies for misf exponential-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# logit specification class membership
cmisfexponormeff_logit <- function(object, level) {
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
  mustar1 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv)) * dnorm(mustar1/sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
  u_c2 <- mustar2 + sqrt(exp(Wv)) * dnorm(mustar2/sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv)) * pnorm(mustar1/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv)) * pnorm(mustar2/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv)) *
      pnorm(mustar1/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv)) *
      pnorm(mustar2/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# cauchit specification class membership
cmisfexponormeff_cauchit <- function(object, level) {
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
  mustar1 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv)) * dnorm(mustar1/sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
  u_c2 <- mustar2 + sqrt(exp(Wv)) * dnorm(mustar2/sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv)) * pnorm(mustar1/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv)) * pnorm(mustar2/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv)) *
      pnorm(mustar1/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv)) *
      pnorm(mustar2/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# probit specification class membership
cmisfexponormeff_probit <- function(object, level) {
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
  mustar1 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv)) * dnorm(mustar1/sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
  u_c2 <- mustar2 + sqrt(exp(Wv)) * dnorm(mustar2/sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv)) * pnorm(mustar1/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv)) * pnorm(mustar2/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv)) *
      pnorm(mustar1/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv)) *
      pnorm(mustar2/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# cloglog specification class membership
cmisfexponormeff_cloglog <- function(object, level) {
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
  mustar1 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv)) * dnorm(mustar1/sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
  u_c2 <- mustar2 + sqrt(exp(Wv)) * dnorm(mustar2/sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv)) * pnorm(mustar1/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv)) * pnorm(mustar2/sqrt(exp(Wv)) -
      sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv)) *
      pnorm(mustar1/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar1/sqrt(exp(Wv)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv)) *
      pnorm(mustar2/sqrt(exp(Wv)) + sqrt(exp(Wv)))/pnorm(mustar2/sqrt(exp(Wv)))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for misf exponential-normal distribution
#' @param object object of class sfacross
#' @noRd
# logit specification class membership
cmisfmargexponorm_Eu_logit <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmargexponorm_Vu_logit <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# cauchit specification class membership
cmisfmargexponorm_Eu_cauchit <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmargexponorm_Vu_cauchit <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# probit specification class membership
cmisfmargexponorm_Eu_probit <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmargexponorm_Vu_probit <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# cloglog specification class membership
cmisfmargexponorm_Eu_cloglog <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmisfmargexponorm_Vu_cloglog <- function(object) {
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
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv/2) -
    exp(Wv/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}
