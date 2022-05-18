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
# Convolution: exponential - normal                                            #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf exponential-normal distribution
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
# same sigma_u

## logit specification class membership
ccnsfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
ccnsfexponormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
ccnsfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
ccnsfexponormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# different sigma_u

## logit specification class membership
cmcesfexponormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmcesfexponormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmcesfexponormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmcesfexponormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/exp(Wu1/2) * exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf exponential-normal distribution
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
# same sigma_u
cstcnsfexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + exponential - normal distributions...\n")
  initExpo <- maxLik(logLik = cexponormlike, start = cstexponorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradexponormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initExpo$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("CN_", colnames(Zvar)))
  names(initExpo$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initExpo = initExpo))
}

# different sigma_u
cstmcesfexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
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
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("MCE_",
    colnames(Zvar)))
  names(initExpo$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initExpo = initExpo))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf exponential-normal distribution
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
# same sigma_u

## logit specification class membership
cgradcnsfexponormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
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
  musig1 <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu_h)
  sigx2_1 <- exp(exp(Wv1)/(2 * exp(Wu)) + S * (epsilon)/ewu_h)
  sigx2_2 <- exp(exp(Wv2)/(2 * exp(Wu)) + S * (epsilon)/ewu_h)
  sigx3 <- (prC * sigx1_2 * sigx2_2 + sigx1_1 * sigx2_1 * ewz/wzdeno)
  sigx4 <- (prC * sigx2_2 * pmusig2 + sigx2_1 * ewz * pmusig1/wzdeno)
  sigx5_1 <- (dmusig1 * ewv1_h/ewu_h)
  sigx5_2 <- (dmusig2 * ewv2_h/ewu_h)
  epsiu <- (S * (epsilon)/ewu_h)
  sigx6_1 <- (exp(Wu) * exp(Wv1)/(2 * exp(Wu))^2)
  sigx6_2 <- (exp(Wu) * exp(Wv2)/(2 * exp(Wu))^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (sigx7_1 * sigx2_1 * ewz/wzdeno + sigx7_2 * prC *
    sigx2_2)
  epsiv1 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  epsiv2 <- (0.5 * (ewv2_h/ewu_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx9_1 <- (exp(Wv1) * pmusig1/(2 * exp(Wu)) - epsiv1 * dmusig1)
  sigx9_2 <- (exp(Wv2) * pmusig2/(2 * exp(Wu)) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  swz <- (sigx4 * wzdeno)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx8/sigx4 -
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx2_1 *
    sigx9_1 * ewz/swz, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = prC * sigx2_2 * sigx9_2/sigx4, FUN = "*"), sweep(Zvar,
    MARGIN = 1, STATS = prC * sigx10 * ewz/swz, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradcnsfexponormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3 <- (ewz2 * sigx1_2 * sigx2_2 + ewz1 * sigx1_1 * sigx2_1)
  sigx4 <- (ewz2 * sigx2_2 * pmusig2 + ewz1 * sigx2_1 * pmusig1)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr)
  epsiu <- (S * (epsilon)/ewusr)
  sigx6_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx6_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (ewz2 * sigx7_2 * sigx2_2 + sigx7_1 * ewz1 * sigx2_1)
  epsiv1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx8/sigx4 -
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz1 *
    sigx2_1 * sigx9_1/sigx4, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = ewz2 * sigx2_2 * sigx9_2/sigx4, FUN = "*"), sweep(Zvar,
    MARGIN = 1, STATS = sigx10/(pi * sigx4 * ((Wz)^2 + 1)),
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradcnsfexponormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3 <- ((1 - pwZ) * sigx1_2 * sigx2_2 + sigx1_1 * sigx2_1 *
    pwZ)
  sigx4 <- ((1 - pwZ) * sigx2_2 * pmusig2 + sigx2_1 * pmusig1 *
    pwZ)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr)
  epsiu <- (S * (epsilon)/ewusr)
  sigx6_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx6_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (sigx7_1 * sigx2_1 * pwZ + sigx7_2 * (1 - pwZ) *
    sigx2_2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx8/sigx4 -
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx2_1 *
    sigx9_1 * pwZ/sigx4, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = (1 - pwZ) * sigx2_2 * sigx9_2/sigx4, FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = dwZ * sigx10/sigx4, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradcnsfexponormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3 <- ((1 - prZ) * sigx1_1 * sigx2_1 + sigx1_2 * prZ *
    sigx2_2)
  sigx4 <- ((1 - prZ) * sigx2_1 * pmusig1 + prZ * sigx2_2 *
    pmusig2)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr)
  epsiu <- (S * (epsilon)/ewusr)
  sigx6_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx6_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (sigx7_1 * (1 - prZ) * sigx2_1 + sigx7_2 * prZ *
    sigx2_2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx8/sigx4 -
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (1 -
    prZ) * sigx2_1 * sigx9_1/sigx4, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = prZ * sigx2_2 * sigx9_2/sigx4, FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = prZ * sigx10 * ewz/sigx4,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u

## logit specification class membership
cgradmcesfexponormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
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
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx1_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx2_1 <- (dmusig1/ewv1_h - pmusig1/ewu1_h)
  sigx2_2 <- (dmusig2/ewv2_h - pmusig2/ewu2_h)
  wzu1 <- (wzdeno * ewu1_h)
  wzu2 <- (wzdeno * ewu2_h)
  sigx3 <- (prC * sigx2_2 * sigx1_2/ewu2_h + sigx2_1 * sigx1_1 *
    ewz/wzu1)
  sigx4 <- (prC * sigx1_2 * pmusig2/ewu2_h + sigx1_1 * ewz *
    pmusig1/wzu1)
  duv1 <- (dmusig1 * ewv1_h/ewu1_h)
  duv2 <- (dmusig2 * ewv2_h/ewu2_h)
  ueps1 <- (S * (epsilon)/ewu1_h)
  ueps2 <- (S * (epsilon)/ewu2_h)
  usq1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  usq2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx5_1 <- (0.5 * duv1 - (0.5 * ueps1 + 2 * usq1) * pmusig1)
  sigx5_2 <- (0.5 * duv2 - (0.5 + 0.5 * ueps2 + 2 * usq2) *
    pmusig2)
  sigx6_1 <- (wzdeno * ewu1_h * pmusig1/wzu1^2)
  sigx7_1 <- (sigx5_1/wzu1 - 0.5 * sigx6_1)
  vu1 <- (ewv1_h/ewu1_h)
  vu2 <- (ewv2_h/ewu2_h)
  epsiv1 <- (S * (epsilon)/ewv1_h)
  epsiv2 <- (S * (epsilon)/ewv2_h)
  sigx8_1 <- (ewv1 * pmusig1/(2 * ewu1) - (0.5 * vu1 - 0.5 *
    epsiv1) * dmusig1)
  sigx8_2 <- (ewv2 * pmusig2/(2 * ewu2) - (0.5 * vu2 - 0.5 *
    epsiv2) * dmusig2)
  wzsq <- (1/wzu1 - ewu1_h * ewz/wzu1^2)
  sigx9 <- (wzsq * sigx1_1 * pmusig1 - prC * sigx1_2 * pmusig2/wzu2)
  s4u1 <- (sigx4 * wzdeno * ewu1_h)
  s4u2 <- (sigx4 * ewu2_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx7_1 *
    sigx1_1 * ewz/sigx4, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = sigx5_2 * prC * sigx1_2/s4u2, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = sigx1_1 * sigx8_1 * ewz/s4u1, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = prC * sigx1_2 * sigx8_2/s4u2,
      FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9 *
      ewz/sigx4, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmcesfexponormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- (ewz2 * sigx1_2 * sigx2_2/ewusr2 + ewz1 * sigx1_1 *
    sigx2_1/ewusr1)
  sigx4 <- (ewz2 * sigx2_2 * pmusig2/ewusr2 + ewz1 * sigx2_1 *
    pmusig1/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx8_1 *
    ewz1 * sigx2_1/(sigx4 * ewusr1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = ewz2 * sigx8_2 * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz1 *
    sigx2_1 * sigx9_1/(sigx4 * ewusr1), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = ewz2 * sigx2_2 * sigx9_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx10/(pi *
    sigx4 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmcesfexponormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - pwZ) * sigx1_2 * sigx2_2/ewusr2 + sigx1_1 *
    sigx2_1 * pwZ/ewusr1)
  sigx4 <- ((1 - pwZ) * sigx2_2 * pmusig2/ewusr2 + sigx2_1 *
    pmusig1 * pwZ/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx8_1 *
    pwZ * sigx2_1/(sigx4 * ewusr1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = sigx8_2 * (1 - pwZ) * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx2_1 *
    sigx9_1 * pwZ/(sigx4 * ewusr1), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = (1 - pwZ) * sigx2_2 * sigx9_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = dwZ *
    sigx10/sigx4, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmcesfexponormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - prZ) * sigx1_1 * sigx2_1/ewusr1 + sigx1_2 *
    prZ * sigx2_2/ewusr2)
  sigx4 <- ((1 - prZ) * sigx2_1 * pmusig1/ewusr1 + prZ * sigx2_2 *
    pmusig2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx8_1 *
    (1 - prZ) * sigx2_1/(sigx4 * ewusr1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = prZ * sigx8_2 * sigx2_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (1 -
    prZ) * sigx2_1 * sigx9_1/(sigx4 * ewusr1), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = prZ * sigx2_2 * sigx9_2/(sigx4 *
      ewusr2), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = prZ *
      sigx10 * ewz/sigx4, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf exponential-normal distribution
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
# same sigma_u

## logit specification class membership
chesscnsfexponormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
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
  musig1 <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu_h)
  sigx2_1 <- exp(exp(Wv1)/(2 * exp(Wu)) + S * (epsilon)/ewu_h)
  sigx2_2 <- exp(exp(Wv2)/(2 * exp(Wu)) + S * (epsilon)/ewu_h)
  sigx3 <- (prC * sigx1_2 * sigx2_2 + sigx1_1 * sigx2_1 * ewz/wzdeno)
  sigx4 <- (prC * sigx2_2 * pmusig2 + sigx2_1 * ewz * pmusig1/wzdeno)
  sigx5_1 <- (dmusig1 * ewv1_h/ewu_h)
  sigx5_2 <- (dmusig2 * ewv2_h/ewu_h)
  epsiu <- (S * (epsilon)/ewu_h)
  sigx6_1 <- (exp(Wu) * exp(Wv1)/(2 * exp(Wu))^2)
  sigx6_2 <- (exp(Wu) * exp(Wv2)/(2 * exp(Wu))^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (sigx7_1 * sigx2_1 * ewz/wzdeno + sigx7_2 * prC *
    sigx2_2)
  epsiv1 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  epsiv2 <- (0.5 * (ewv2_h/ewu_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx9_1 <- (exp(Wv1) * pmusig1/(2 * exp(Wu)) - epsiv1 * dmusig1)
  sigx9_2 <- (exp(Wv2) * pmusig2/(2 * exp(Wu)) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  swz <- (sigx4 * wzdeno)
  wuwusq <- (1 - 8 * (exp(Wu)^2/(2 * exp(Wu))^2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (((musig1/ewv1_h - 1/ewu_h) * dmusig1/ewv1_h -
      sigx1_1/ewu_h) * sigx2_1 * ewz/wzdeno + ((musig2/ewv2_h -
      1/ewu_h) * dmusig2/ewv2_h - sigx1_2/ewu_h) * prC *
      sigx2_2 - sigx3^2/sigx4)/sigx4, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 + 0.5 * epsiu +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewu_h + (0.5 *
      (musig1/ewu_h) - (0.5 * epsiu + 2 * sigx6_1)/ewv1_h) *
      dmusig1) * sigx2_1 * ewz/wzdeno + (((0.5 + 0.5 *
      epsiu + 2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewu_h +
      (0.5 * (musig2/ewu_h) - (0.5 * epsiu + 2 * sigx6_2)/ewv2_h) *
        dmusig2) * prC * sigx2_2 - sigx8 * sigx3/sigx4)/sigx4,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((dmusig1 * (exp(Wv1)/(2 * exp(Wu)) - (epsiv1 * musig1 +
    0.5))/ewv1_h - sigx9_1/ewu_h)/swz - sigx3 * wzdeno *
    sigx9_1/swz^2) * sigx2_1 * ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prC * (dmusig2 * (exp(Wv2)/(2 * exp(Wu)) -
      (epsiv2 * musig2 + 0.5))/ewv2_h - (sigx3/sigx4 +
      1/ewu_h) * sigx9_2) * sigx2_2/sigx4, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx1_1 * sigx2_1 -
      sigx1_2 * sigx2_2)/swz - sigx3 * wzdeno * sigx10/swz^2) *
      prC * ewz, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewv1_h * musig1/ewu_h) - 0.5) - 0.5 *
      (0.5 * epsiu + 2 * sigx6_1)) * dmusig1 * ewv1_h/ewu_h -
      (sigx7_1 * (0.5 * epsiu + 2 * sigx6_1) + (2 * (wuwusq *
        exp(Wu) * exp(Wv1)/(2 * exp(Wu))^2) - 0.25 *
        epsiu) * pmusig1)) * sigx2_1 * ewz/wzdeno + ((0.5 *
      (0.5 * (ewv2_h * musig2/ewu_h) - 0.5) - 0.5 * (0.5 *
      epsiu + 2 * sigx6_2)) * dmusig2 * ewv2_h/ewu_h -
      (sigx7_2 * (0.5 * epsiu + 2 * sigx6_2) + (2 * (wuwusq *
        exp(Wu) * exp(Wv2)/(2 * exp(Wu))^2) - 0.25 *
        epsiu) * pmusig2)) * prC * sigx2_2 - sigx8^2/sigx4)/sigx4,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((dmusig1 * ewv1_h/(4 *
      (exp(Wu) * ewu_h)) - 2 * (exp(Wu) * pmusig1/(2 *
      exp(Wu))^2)) * exp(Wv1) - ((0.5 * (epsiv1 * musig1) -
      0.25) * dmusig1 * ewv1_h/ewu_h + (0.5 * epsiu + 2 *
      sigx6_1) * sigx9_1))/swz - sigx8 * wzdeno * sigx9_1/swz^2) *
      sigx2_1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewv2_h/(4 * (exp(Wu) *
      ewu_h)) - 2 * (exp(Wu) * pmusig2/(2 * exp(Wu))^2)) *
      exp(Wv2) - ((sigx8/sigx4 + 0.5 * epsiu + 2 * sigx6_2) *
      sigx9_2 + (0.5 * (epsiv2 * musig2) - 0.25) * dmusig2 *
      ewv2_h/ewu_h)) * prC * sigx2_2/sigx4, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sigx7_1 * sigx2_1 - sigx7_2 *
      sigx2_2)/swz - sigx8 * wzdeno * sigx10/swz^2) * prC *
      ewz, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * exp(Wv1)/exp(Wu) - (0.25 *
      (ewv1_h/ewu_h) + 0.25 * (S * (epsilon)/ewv1_h) -
      epsiv1^2 * musig1) * dmusig1)/swz - sigx2_1 * sigx9_1^2 *
      ewz/swz^2) * sigx2_1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (prC * sigx2_1 * sigx2_2 * sigx9_1 * sigx9_2 * ewz/(sigx4^2 *
      wzdeno)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (1/swz - sigx10 * ewz/swz^2) *
      sigx2_1 * sigx9_1 * ewz, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx9_2/2 + (pmusig2 - epsiv2 * dmusig2)/2) *
      exp(Wv2)/exp(Wu) - ((0.25 * (ewv2_h/ewu_h) + 0.25 *
      (S * (epsilon)/ewv2_h) - epsiv2^2 * musig2) * dmusig2 +
      prC * sigx2_2 * sigx9_2^2/sigx4)) * prC * sigx2_2/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((prC * wzdeno * sigx10/swz^2 +
      1/swz) * prC * sigx2_2 * sigx9_2 * ewz), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (prC/swz - sigx2_1 * ewz *
      pmusig1/swz^2) * prC * sigx10 * ewz, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesscnsfexponormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3 <- (ewz2 * sigx1_2 * sigx2_2 + ewz1 * sigx1_1 * sigx2_1)
  sigx4 <- (ewz2 * sigx2_2 * pmusig2 + ewz1 * sigx2_1 * pmusig1)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr)
  epsiu <- (S * (epsilon)/ewusr)
  sigx6_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx6_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (ewz2 * sigx7_2 * sigx2_2 + sigx7_1 * ewz1 * sigx2_1)
  epsiv1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr1 - 1/ewusr) * dmusig1/ewvsr1 -
      sigx1_1/ewusr) * ewz1 * sigx2_1 + ((musig2/ewvsr2 -
      1/ewusr) * dmusig2/ewvsr2 - sigx1_2/ewusr) * ewz2 *
      sigx2_2 - sigx3^2/sigx4)/sigx4, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 + 0.5 * epsiu +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr + (0.5 *
      (musig1/ewusr) - (0.5 * epsiu + 2 * sigx6_1)/ewvsr1) *
      dmusig1) * ewz1 * sigx2_1 + (((0.5 + 0.5 * epsiu +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr + (0.5 *
      (musig2/ewusr) - (0.5 * epsiu + 2 * sigx6_2)/ewvsr2) *
      dmusig2) * ewz2 * sigx2_2 - sigx8 * sigx3/sigx4)/sigx4,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ewz1 * (dmusig1 * (ewv1/(2 * ewu) - (epsiv1 * musig1 +
    0.5))/ewvsr1 - (sigx3/sigx4 + 1/ewusr) * sigx9_1) * sigx2_1/sigx4,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ewz2 * (dmusig2 * (ewv2/(2 * ewu) -
      (epsiv2 * musig2 + 0.5))/ewvsr2 - (sigx3/sigx4 +
      1/ewusr) * sigx9_2) * sigx2_2/sigx4, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx1_1 * sigx2_1 -
      sigx1_2 * sigx2_2)/(pi * sigx4 * ((Wz)^2 + 1)) -
      pi * sigx3 * ((Wz)^2 + 1) * sigx10/(pi * sigx4 *
        ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr) - 0.5) - 0.5 *
      (0.5 * epsiu + 2 * sigx6_1)) * dmusig1 * ewvsr1/ewusr -
      (sigx7_1 * (0.5 * epsiu + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) -
        0.25 * epsiu) * pmusig1)) * ewz1 * sigx2_1 +
      ((0.5 * (0.5 * (ewvsr2 * musig2/ewusr) - 0.5) - 0.5 *
        (0.5 * epsiu + 2 * sigx6_2)) * dmusig2 * ewvsr2/ewusr -
        (sigx7_2 * (0.5 * epsiu + 2 * sigx6_2) + (2 *
          ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv2/(2 *
          ewu)^2) - 0.25 * epsiu) * pmusig2)) * ewz2 *
        sigx2_2 - sigx8^2/sigx4)/sigx4, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig1 * ewvsr1/(4 * (ewu *
      ewusr)) - 2 * (ewu * pmusig1/(2 * ewu)^2)) * ewv1 -
      ((sigx8/sigx4 + 0.5 * epsiu + 2 * sigx6_1) * sigx9_1 +
        (0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 *
          ewvsr1/ewusr)) * ewz1 * sigx2_1/sigx4, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewvsr2/(4 * (ewu *
      ewusr)) - 2 * (ewu * pmusig2/(2 * ewu)^2)) * ewv2 -
      ((sigx8/sigx4 + 0.5 * epsiu + 2 * sigx6_2) * sigx9_2 +
        (0.5 * (epsiv2 * musig2) - 0.25) * dmusig2 *
          ewvsr2/ewusr)) * ewz2 * sigx2_2/sigx4, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sigx7_1 * sigx2_1 - sigx7_2 *
      sigx2_2)/(pi * sigx4 * ((Wz)^2 + 1)) - pi * sigx8 *
      ((Wz)^2 + 1) * sigx10/(pi * sigx4 * ((Wz)^2 + 1))^2),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv1/ewu - ((0.25 * (ewvsr1/ewusr) +
      0.25 * (S * (epsilon)/ewvsr1) - epsiv1^2 * musig1) *
      dmusig1 + ewz1 * sigx2_1 * sigx9_1^2/sigx4)) * ewz1 *
      sigx2_1/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (ewz2 * ewz1 * sigx2_1 * sigx2_2 * sigx9_1 * sigx9_2/sigx4^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1/(pi * sigx4 * ((Wz)^2 +
      1)) - pi * ((Wz)^2 + 1) * ewz1 * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2) * sigx2_1 * sigx9_1, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx9_2/2 + (pmusig2 - epsiv2 * dmusig2)/2) *
      ewv2/ewu - ((0.25 * (ewvsr2/ewusr) + 0.25 * (S *
      (epsilon)/ewvsr2) - epsiv2^2 * musig2) * dmusig2 +
      ewz2 * sigx2_2 * sigx9_2^2/sigx4)) * ewz2 * sigx2_2/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1/(pi * sigx4 * ((Wz)^2 +
      1)) + pi * ((Wz)^2 + 1) * ewz2 * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2) * sigx2_2 * sigx9_2), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx4) +
      sigx2_1 * pmusig1 - sigx2_2 * pmusig2) * sigx10/(pi *
      sigx4 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesscnsfexponormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3 <- ((1 - pwZ) * sigx1_2 * sigx2_2 + sigx1_1 * sigx2_1 *
    pwZ)
  sigx4 <- ((1 - pwZ) * sigx2_2 * pmusig2 + sigx2_1 * pmusig1 *
    pwZ)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr)
  epsiu <- (S * (epsilon)/ewusr)
  sigx6_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx6_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (sigx7_1 * sigx2_1 * pwZ + sigx7_2 * (1 - pwZ) *
    sigx2_2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr1 - 1/ewusr) * dmusig1/ewvsr1 -
      sigx1_1/ewusr) * sigx2_1 * pwZ + ((musig2/ewvsr2 -
      1/ewusr) * dmusig2/ewvsr2 - sigx1_2/ewusr) * (1 -
      pwZ) * sigx2_2 - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 + 0.5 * epsiu +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr + (0.5 *
      (musig1/ewusr) - (0.5 * epsiu + 2 * sigx6_1)/ewvsr1) *
      dmusig1) * sigx2_1 * pwZ + (((0.5 + 0.5 * epsiu +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr + (0.5 *
      (musig2/ewusr) - (0.5 * epsiu + 2 * sigx6_2)/ewvsr2) *
      dmusig2) * (1 - pwZ) * sigx2_2 - sigx8 * sigx3/sigx4)/sigx4,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (dmusig1 * (ewv1/(2 * ewu) - (epsiv1 * musig1 + 0.5))/ewvsr1 -
    (sigx3/sigx4 + 1/ewusr) * sigx9_1) * sigx2_1 * pwZ/sigx4,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (1 - pwZ) * (dmusig2 * (ewv2/(2 *
      ewu) - (epsiv2 * musig2 + 0.5))/ewvsr2 - (sigx3/sigx4 +
      1/ewusr) * sigx9_2) * sigx2_2/sigx4, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1 * sigx2_1 -
      (sigx3 * sigx10/sigx4 + sigx1_2 * sigx2_2)) * dwZ/sigx4,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr) - 0.5) - 0.5 *
      (0.5 * epsiu + 2 * sigx6_1)) * dmusig1 * ewvsr1/ewusr -
      (sigx7_1 * (0.5 * epsiu + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) -
        0.25 * epsiu) * pmusig1)) * sigx2_1 * pwZ + ((0.5 *
      (0.5 * (ewvsr2 * musig2/ewusr) - 0.5) - 0.5 * (0.5 *
      epsiu + 2 * sigx6_2)) * dmusig2 * ewvsr2/ewusr -
      (sigx7_2 * (0.5 * epsiu + 2 * sigx6_2) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv2/(2 * ewu)^2) -
        0.25 * epsiu) * pmusig2)) * (1 - pwZ) * sigx2_2 -
      sigx8^2/sigx4)/sigx4, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig1 * ewvsr1/(4 * (ewu *
      ewusr)) - 2 * (ewu * pmusig1/(2 * ewu)^2)) * ewv1 -
      ((sigx8/sigx4 + 0.5 * epsiu + 2 * sigx6_1) * sigx9_1 +
        (0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 *
          ewvsr1/ewusr)) * sigx2_1 * pwZ/sigx4, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewvsr2/(4 * (ewu *
      ewusr)) - 2 * (ewu * pmusig2/(2 * ewu)^2)) * ewv2 -
      ((sigx8/sigx4 + 0.5 * epsiu + 2 * sigx6_2) * sigx9_2 +
        (0.5 * (epsiv2 * musig2) - 0.25) * dmusig2 *
          ewvsr2/ewusr)) * (1 - pwZ) * sigx2_2/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx7_1 * sigx2_1 - (sigx8 *
      sigx10/sigx4 + sigx7_2 * sigx2_2)) * dwZ/sigx4, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv1/ewu - ((0.25 * (ewvsr1/ewusr) +
      0.25 * (S * (epsilon)/ewvsr1) - epsiv1^2 * musig1) *
      dmusig1 + sigx2_1 * sigx9_1^2 * pwZ/sigx4)) * sigx2_1 *
      pwZ/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((1 - pwZ) * sigx2_1 * sigx2_2 * sigx9_1 * sigx9_2 *
      pwZ/sigx4^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx10 * pwZ/sigx4) *
      dwZ * sigx2_1 * sigx9_1/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx9_2/2 + (pmusig2 - epsiv2 * dmusig2)/2) *
      ewv2/ewu - ((0.25 * (ewvsr2/ewusr) + 0.25 * (S *
      (epsilon)/ewvsr2) - epsiv2^2 * musig2) * dmusig2 +
      (1 - pwZ) * sigx2_2 * sigx9_2^2/sigx4)) * (1 - pwZ) *
      sigx2_2/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((1 - pwZ) * sigx10/sigx4 +
      1) * dwZ * sigx2_2 * sigx9_2/sigx4), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * (dwZ * (dwZ * sigx10/sigx4 +
      Wz) * sigx10/sigx4), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesscnsfexponormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewu <- exp(Wu)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3 <- ((1 - prZ) * sigx1_1 * sigx2_1 + sigx1_2 * prZ *
    sigx2_2)
  sigx4 <- ((1 - prZ) * sigx2_1 * pmusig1 + prZ * sigx2_2 *
    pmusig2)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr)
  epsiu <- (S * (epsilon)/ewusr)
  sigx6_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx6_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu + 2 * sigx6_2) *
    pmusig2)
  sigx8 <- (sigx7_1 * (1 - prZ) * sigx2_1 + sigx7_2 * prZ *
    sigx2_2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1 - sigx2_2 * pmusig2)
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr1 - 1/ewusr) * dmusig1/ewvsr1 -
      sigx1_1/ewusr) * (1 - prZ) * sigx2_1 + ((musig2/ewvsr2 -
      1/ewusr) * dmusig2/ewvsr2 - sigx1_2/ewusr) * prZ *
      sigx2_2 - sigx3^2/sigx4)/sigx4, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 + 0.5 * epsiu +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr + (0.5 *
      (musig1/ewusr) - (0.5 * epsiu + 2 * sigx6_1)/ewvsr1) *
      dmusig1) * (1 - prZ) * sigx2_1 + (((0.5 + 0.5 * epsiu +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr + (0.5 *
      (musig2/ewusr) - (0.5 * epsiu + 2 * sigx6_2)/ewvsr2) *
      dmusig2) * prZ * sigx2_2 - sigx8 * sigx3/sigx4)/sigx4,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (1 - prZ) * (dmusig1 * (ewv1/(2 * ewu) - (epsiv1 *
    musig1 + 0.5))/ewvsr1 - (sigx3/sigx4 + 1/ewusr) * sigx9_1) *
    sigx2_1/sigx4, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (dmusig2 * (ewv2/(2 * ewu) - (epsiv2 *
      musig2 + 0.5))/ewvsr2 - (sigx3/sigx4 + 1/ewusr) *
      sigx9_2) * prZ * sigx2_2/sigx4, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1 * sigx2_1 -
      (sigx3 * sigx10/sigx4 + sigx1_2 * sigx2_2)) * prZ *
      ewz/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr) - 0.5) - 0.5 *
      (0.5 * epsiu + 2 * sigx6_1)) * dmusig1 * ewvsr1/ewusr -
      (sigx7_1 * (0.5 * epsiu + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) -
        0.25 * epsiu) * pmusig1)) * (1 - prZ) * sigx2_1 +
      ((0.5 * (0.5 * (ewvsr2 * musig2/ewusr) - 0.5) - 0.5 *
        (0.5 * epsiu + 2 * sigx6_2)) * dmusig2 * ewvsr2/ewusr -
        (sigx7_2 * (0.5 * epsiu + 2 * sigx6_2) + (2 *
          ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv2/(2 *
          ewu)^2) - 0.25 * epsiu) * pmusig2)) * prZ *
        sigx2_2 - sigx8^2/sigx4)/sigx4, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig1 * ewvsr1/(4 * (ewu *
      ewusr)) - 2 * (ewu * pmusig1/(2 * ewu)^2)) * ewv1 -
      ((sigx8/sigx4 + 0.5 * epsiu + 2 * sigx6_1) * sigx9_1 +
        (0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 *
          ewvsr1/ewusr)) * (1 - prZ) * sigx2_1/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((dmusig2 * ewvsr2/(4 * (ewu *
      ewusr)) - 2 * (ewu * pmusig2/(2 * ewu)^2)) * ewv2 -
      ((sigx8/sigx4 + 0.5 * epsiu + 2 * sigx6_2) * sigx9_2 +
        (0.5 * (epsiv2 * musig2) - 0.25) * dmusig2 *
          ewvsr2/ewusr)) * prZ * sigx2_2/sigx4, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx7_1 * sigx2_1 - (sigx8 *
      sigx10/sigx4 + sigx7_2 * sigx2_2)) * prZ * ewz/sigx4,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv1/ewu - ((0.25 * (ewvsr1/ewusr) +
      0.25 * (S * (epsilon)/ewvsr1) - epsiv1^2 * musig1) *
      dmusig1 + (1 - prZ) * sigx2_1 * sigx9_1^2/sigx4)) *
      (1 - prZ) * sigx2_1/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((1 - prZ) * prZ * sigx2_1 * sigx2_2 * sigx9_1 * sigx9_2/sigx4^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 - prZ) * sigx10/sigx4) *
      prZ * sigx2_1 * sigx9_1 * ewz/sigx4, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx9_2/2 + (pmusig2 - epsiv2 * dmusig2)/2) *
      ewv2/ewu - ((0.25 * (ewvsr2/ewusr) + 0.25 * (S *
      (epsilon)/ewvsr2) - epsiv2^2 * musig2) * dmusig2 +
      prZ * sigx2_2 * sigx9_2^2/sigx4)) * prZ * sigx2_2/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1 + prZ * sigx10/sigx4) *
      prZ * sigx2_2 * sigx9_2 * ewz/sigx4), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 + prZ * sigx10/sigx4) *
      ewz) * prZ * sigx10 * ewz/sigx4, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# different sigma_u

## logit specification class membership
chessmcesfexponormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
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
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx1_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx2_1 <- (dmusig1/ewv1_h - pmusig1/ewu1_h)
  sigx2_2 <- (dmusig2/ewv2_h - pmusig2/ewu2_h)
  wzu1 <- (wzdeno * ewu1_h)
  wzu2 <- (wzdeno * ewu2_h)
  sigx3 <- (prC * sigx2_2 * sigx1_2/ewu2_h + sigx2_1 * sigx1_1 *
    ewz/wzu1)
  sigx4 <- (prC * sigx1_2 * pmusig2/ewu2_h + sigx1_1 * ewz *
    pmusig1/wzu1)
  duv1 <- (dmusig1 * ewv1_h/ewu1_h)
  duv2 <- (dmusig2 * ewv2_h/ewu2_h)
  ueps1 <- (S * (epsilon)/ewu1_h)
  ueps2 <- (S * (epsilon)/ewu2_h)
  usq1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  usq2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx5_1 <- (0.5 * duv1 - (0.5 * ueps1 + 2 * usq1) * pmusig1)
  sigx5_2 <- (0.5 * duv2 - (0.5 + 0.5 * ueps2 + 2 * usq2) *
    pmusig2)
  sigx6_1 <- (wzdeno * ewu1_h * pmusig1/wzu1^2)
  sigx7_1 <- (sigx5_1/wzu1 - 0.5 * sigx6_1)
  vu1 <- (ewv1_h/ewu1_h)
  vu2 <- (ewv2_h/ewu2_h)
  epsiv1 <- (S * (epsilon)/ewv1_h)
  epsiv2 <- (S * (epsilon)/ewv2_h)
  sigx8_1 <- (ewv1 * pmusig1/(2 * ewu1) - (0.5 * vu1 - 0.5 *
    epsiv1) * dmusig1)
  sigx8_2 <- (ewv2 * pmusig2/(2 * ewu2) - (0.5 * vu2 - 0.5 *
    epsiv2) * dmusig2)
  wzsq <- (1/wzu1 - ewu1_h * ewz/wzu1^2)
  sigx9 <- (wzsq * sigx1_1 * pmusig1 - prC * sigx1_2 * pmusig2/wzu2)
  s4u1 <- (sigx4 * wzdeno * ewu1_h)
  s4u2 <- (sigx4 * ewu2_h)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (((musig1/ewv1_h - 1/ewu1_h) *
      dmusig1/ewv1_h - sigx2_1/ewu1_h) * sigx1_1 * ewz/wzu1 +
      ((musig2/ewv2_h - 1/ewu2_h) * dmusig2/ewv2_h - sigx2_2/ewu2_h) *
        prC * sigx1_2/ewu2_h - sigx3^2/sigx4)/sigx4,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((((0.5 + 0.5 * ueps1 +
      2 * usq1) * pmusig1 - 0.5 * duv1)/ewu1_h + (0.5 *
      (musig1/ewu1_h) - (0.5 * ueps1 + 2 * usq1)/ewv1_h) *
      dmusig1)/wzdeno + 0.5 * sigx6_1)/ewu1_h - (sigx7_1 *
      sigx3/sigx4 + 0.5 * (wzdeno * dmusig1 * ewu1_h/(wzu1^2 *
      ewv1_h)))) * sigx1_1 * ewz/sigx4, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ueps2 + 1 +
      2 * usq2) * pmusig2 - 0.5 * duv2)/ewu2_h + (0.5 *
      (musig2/ewu2_h) - (0.5 + 0.5 * ueps2 + 2 * usq2)/ewv2_h) *
      dmusig2)/s4u2 - sigx3 * sigx5_2 * ewu2_h/s4u2^2) *
      prC * sigx1_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((dmusig1 * (ewv1/(2 * ewu1) - ((0.5 * vu1 - 0.5 *
    epsiv1) * musig1 + 0.5))/ewv1_h - sigx8_1/ewu1_h)/s4u1 -
    sigx3 * wzdeno * ewu1_h * sigx8_1/s4u1^2) * sigx1_1 *
    ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmusig2 * (ewv2/(2 *
      ewu2) - ((0.5 * vu2 - 0.5 * epsiv2) * musig2 + 0.5))/ewv2_h -
      sigx8_2/ewu2_h)/s4u2 - sigx3 * ewu2_h * sigx8_2/s4u2^2) *
      prC * sigx1_2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (wzsq * sigx2_1 * sigx1_1 -
      (sigx3 * sigx9/sigx4 + prC * sigx2_2 * sigx1_2/wzu2)) *
      ewz/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewv1_h * musig1/ewu1_h) - 0.5) - 0.5 *
      (0.5 * ueps1 + 2 * usq1)) * dmusig1 * ewv1_h/ewu1_h -
      (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv1/(2 *
        ewu1)^2) - 0.25 * ueps1) * pmusig1)/wzu1 - ((sigx7_1 *
      sigx1_1 * ewz/sigx4 + 0.5 * ueps1 + 2 * usq1) * sigx7_1 +
      (0.5 * ((0.5 - wzdeno^2 * ewu1_h^2/wzu1^2) * ewu1_h *
        pmusig1 + 0.5 * (dmusig1 * ewv1_h)) + 0.5 * (sigx5_1 *
        ewu1_h)) * wzdeno/wzu1^2)) * sigx1_1 * ewz/sigx4,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx7_1 * sigx5_2 * prC * ewu2_h *
      sigx1_1 * sigx1_2 * ewz/s4u2^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((dmusig1 * ewv1_h/(4 *
      (ewu1 * ewu1_h)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * ((0.5 * vu1 - 0.5 * epsiv1) * musig1) -
      0.25) * dmusig1 * ewv1_h/ewu1_h + (0.5 * ueps1 +
      2 * usq1) * sigx8_1))/s4u1 - (sigx7_1 * sigx1_1 *
      ewz + 0.5 * sigx4) * wzdeno * ewu1_h * sigx8_1/s4u1^2) *
      sigx1_1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx7_1 * prC * ewu2_h *
      sigx1_1 * sigx1_2 * sigx8_2 * ewz/s4u2^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (wzsq * dmusig1 * ewv1_h/ewu1_h) - ((((0.5 - wzdeno^2 *
      ewu1_h^2/wzu1^2) * ewz + 0.5 * wzdeno) * ewu1_h/wzu1^2 +
      (0.5 * ueps1 + 2 * usq1) * wzsq) * pmusig1 + sigx7_1 *
      sigx9 * ewz/sigx4)) * sigx1_1 * ewz/sigx4, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewv2_h *
      musig2/ewu2_h) - 0.5) - 0.5 * (0.5 + 0.5 * ueps2 +
      2 * usq2)) * dmusig2 * ewv2_h/ewu2_h - (sigx5_2 *
      (0.5 * ueps2 + 2 * usq2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 * ueps2) *
      pmusig2))/s4u2 - (sigx5_2 * prC * sigx1_2 + 0.5 *
      s4u2) * sigx5_2/s4u2^2) * prC * sigx1_2, FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx5_2 * prC * wzdeno *
      ewu1_h * sigx1_1 * sigx1_2 * sigx8_1 * ewz/(s4u1^2 *
      ewu2_h)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((dmusig2 * ewv2_h/(4 * (ewu2 * ewu2_h)) - 2 * (ewu2 *
      pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * ((0.5 * vu2 -
      0.5 * epsiv2) * musig2) - 0.25) * dmusig2 * ewv2_h/ewu2_h +
      (0.5 * ueps2 + 2 * usq2) * sigx8_2))/s4u2 - (sigx5_2 *
      prC * sigx1_2 + 0.5 * s4u2) * sigx8_2/s4u2^2) * prC *
    sigx1_2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (((sigx9 * sigx5_2/sigx4 + (0.5 * duv2 -
      (0.5 * ueps2 + 2 * usq2) * pmusig2)/wzdeno)/ewu2_h -
      0.5 * (wzdeno * ewu2_h * pmusig2/wzu2^2)) * prC *
      sigx1_2 * ewz/sigx4), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx8_1/2 + (pmusig1 -
      (0.5 * vu1 - 0.5 * epsiv1) * dmusig1)/2) * ewv1/ewu1 -
      (0.25 * vu1 + 0.25 * epsiv1 - (0.5 * vu1 - 0.5 *
        epsiv1)^2 * musig1) * dmusig1)/s4u1 - sigx1_1 *
      sigx8_1^2 * ewz/s4u1^2) * sigx1_1 * ewz, FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (prC * ewu2_h * sigx1_1 * sigx1_2 *
      sigx8_1 * sigx8_2 * ewz/(s4u2^2 * wzdeno * ewu1_h)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1/wzu1 - (sigx9/s4u1 + ewu1_h/wzu1^2) *
      ewz) * sigx1_1 * sigx8_1 * ewz/sigx4, FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx8_2/2 + (pmusig2 -
      (0.5 * vu2 - 0.5 * epsiv2) * dmusig2)/2) * ewv2/ewu2 -
      (0.25 * vu2 + 0.25 * epsiv2 - (0.5 * vu2 - 0.5 *
        epsiv2)^2 * musig2) * dmusig2)/s4u2 - prC * sigx1_2 *
      sigx8_2^2/s4u2^2) * prC * sigx1_2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx9/sigx4 + 1/wzdeno) *
      prC * sigx1_2 * sigx8_2 * ewz/s4u2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    ((prC * (1/(wzdeno^2 * ewu2_h) + ewu2_h/wzu2^2) * sigx1_2 *
      pmusig2 - (sigx9^2/sigx4 + (2 - 2 * (wzdeno * ewu1_h^2 *
      ewz/wzu1^2)) * ewu1_h * sigx1_1 * pmusig1/wzu1^2)) *
      ewz + wzsq * sigx1_1 * pmusig1 - prC * sigx1_2 *
      pmusig2/wzu2) * ewz/sigx4, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmcesfexponormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- (ewz2 * sigx1_2 * sigx2_2/ewusr2 + ewz1 * sigx1_1 *
    sigx2_1/ewusr1)
  sigx4 <- (ewz2 * sigx2_2 * pmusig2/ewusr2 + ewz1 * sigx2_1 *
    pmusig1/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr1 - 1/ewusr1) * dmusig1/ewvsr1 -
      sigx1_1/ewusr1) * ewz1 * sigx2_1/ewusr1 + ((musig2/ewvsr2 -
      1/ewusr2) * dmusig2/ewvsr2 - sigx1_2/ewusr2) * ewz2 *
      sigx2_2/ewusr2 - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu1 + 1 +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 + 0.5 * epsiu1 + 2 *
        sigx6_1)/ewvsr1) * dmusig1)/(sigx4 * ewusr1) -
      sigx3 * sigx8_1 * ewusr1/(sigx4 * ewusr1)^2) * ewz1 *
      sigx2_1, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu2 + 1 +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 + 0.5 * epsiu2 + 2 *
        sigx6_2)/ewvsr2) * dmusig2)/(sigx4 * ewusr2) -
      sigx3 * sigx8_2 * ewusr2/(sigx4 * ewusr2)^2) * ewz2 *
      sigx2_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((dmusig1 * (ewv1/(2 * ewu1) - (epsiv1 * musig1 +
    0.5))/ewvsr1 - sigx9_1/ewusr1)/(sigx4 * ewusr1) - sigx3 *
    ewusr1 * sigx9_1/(sigx4 * ewusr1)^2) * ewz1 * sigx2_1,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmusig2 * (ewv2/(2 *
      ewu2) - (epsiv2 * musig2 + 0.5))/ewvsr2 - sigx9_2/ewusr2)/(sigx4 *
      ewusr2) - sigx3 * ewusr2 * sigx9_2/(sigx4 * ewusr2)^2) *
      ewz2 * sigx2_2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx1_1 * sigx2_1/ewusr1 -
      sigx1_2 * sigx2_2/ewusr2)/(pi * sigx4 * ((Wz)^2 +
      1)) - pi * sigx3 * ((Wz)^2 + 1) * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr1) - 0.5) - 0.5 *
      (0.5 + 0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr1/ewusr1 -
      (sigx8_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv1/(2 *
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
    MARGIN = 1, STATS = wHvar * (((dmusig1 * ewvsr1/(4 *
      (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 *
      ewvsr1/ewusr1 + (0.5 * epsiu1 + 2 * sigx6_1) * sigx9_1))/(sigx4 *
      ewusr1) - (sigx8_1 * ewz1 * sigx2_1 + 0.5 * (sigx4 *
      ewusr1)) * sigx9_1/(sigx4 * ewusr1)^2) * ewz1 * sigx2_1,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * sigx8_1 * ewz1 *
      ewusr2 * sigx2_1 * sigx2_2 * sigx9_2/((sigx4 * ewusr2)^2 *
      ewusr1)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx8_1 * (1/(pi * sigx4 * ((Wz)^2 + 1)) - pi * ((Wz)^2 +
    1) * ewz1 * sigx10/(pi * sigx4 * ((Wz)^2 + 1))^2) * sigx2_1/ewusr1,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr2 *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 + 0.5 * epsiu2 +
      2 * sigx6_2)) * dmusig2 * ewvsr2/ewusr2 - (sigx8_2 *
      (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 * epsiu2) *
      pmusig2))/(sigx4 * ewusr2) - (ewz2 * sigx8_2 * sigx2_2 +
      0.5 * (sigx4 * ewusr2)) * sigx8_2/(sigx4 * ewusr2)^2) *
      ewz2 * sigx2_2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * sigx8_2 * ewz1 *
      ewusr1 * sigx2_1 * sigx2_2 * sigx9_1/((sigx4 * ewusr1)^2 *
      ewusr2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((dmusig2 * ewvsr2/(4 * (ewu2 * ewusr2)) - 2 * (ewu2 *
      pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (epsiv2 *
      musig2) - 0.25) * dmusig2 * ewvsr2/ewusr2 + (0.5 *
      epsiu2 + 2 * sigx6_2) * sigx9_2))/(sigx4 * ewusr2) -
      (ewz2 * sigx8_2 * sigx2_2 + 0.5 * (sigx4 * ewusr2)) *
        sigx9_2/(sigx4 * ewusr2)^2) * ewz2 * sigx2_2,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx8_2 * (1/(pi * sigx4 * ((Wz)^2 +
      1)) + pi * ((Wz)^2 + 1) * ewz2 * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2) * sigx2_2/ewusr2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewvsr1/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr1) - epsiv1^2 * musig1) *
      dmusig1)/(sigx4 * ewusr1) - ewz1 * sigx2_1 * sigx9_1^2/(sigx4 *
      ewusr1)^2) * ewz1 * sigx2_1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (ewz2 * ewz1 * ewusr2 * sigx2_1 * sigx2_2 *
      sigx9_1 * sigx9_2/((sigx4 * ewusr2)^2 * ewusr1)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1/(pi * sigx4 * ((Wz)^2 +
      1)) - pi * ((Wz)^2 + 1) * ewz1 * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2) * sigx2_1 * sigx9_1/ewusr1, FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_2/2 + (pmusig2 -
      epsiv2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewvsr2/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr2) - epsiv2^2 * musig2) *
      dmusig2)/(sigx4 * ewusr2) - ewz2 * sigx2_2 * sigx9_2^2/(sigx4 *
      ewusr2)^2) * ewz2 * sigx2_2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1/(pi * sigx4 * ((Wz)^2 +
      1)) + pi * ((Wz)^2 + 1) * ewz2 * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2) * sigx2_2 * sigx9_2/ewusr2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar *
    ((2 * (pi * Wz * sigx4) + sigx2_1 * pmusig1/ewusr1 -
      sigx2_2 * pmusig2/ewusr2) * sigx10/(pi * sigx4 *
      ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmcesfexponormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - pwZ) * sigx1_2 * sigx2_2/ewusr2 + sigx1_1 *
    sigx2_1 * pwZ/ewusr1)
  sigx4 <- ((1 - pwZ) * sigx2_2 * pmusig2/ewusr2 + sigx2_1 *
    pmusig1 * pwZ/ewusr1)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr1 - 1/ewusr1) * dmusig1/ewvsr1 -
      sigx1_1/ewusr1) * pwZ * sigx2_1/ewusr1 + ((musig2/ewvsr2 -
      1/ewusr2) * dmusig2/ewvsr2 - sigx1_2/ewusr2) * (1 -
      pwZ) * sigx2_2/ewusr2 - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu1 + 1 +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 + 0.5 * epsiu1 + 2 *
        sigx6_1)/ewvsr1) * dmusig1)/(sigx4 * ewusr1) -
      sigx3 * sigx8_1 * ewusr1/(sigx4 * ewusr1)^2) * pwZ *
      sigx2_1, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu2 + 1 +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 + 0.5 * epsiu2 + 2 *
        sigx6_2)/ewvsr2) * dmusig2)/(sigx4 * ewusr2) -
      sigx3 * sigx8_2 * ewusr2/(sigx4 * ewusr2)^2) * (1 -
      pwZ) * sigx2_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((dmusig1 * (ewv1/(2 * ewu1) - (epsiv1 * musig1 +
    0.5))/ewvsr1 - sigx9_1/ewusr1)/(sigx4 * ewusr1) - sigx3 *
    ewusr1 * sigx9_1/(sigx4 * ewusr1)^2) * pwZ * sigx2_1,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmusig2 * (ewv2/(2 *
      ewu2) - (epsiv2 * musig2 + 0.5))/ewvsr2 - sigx9_2/ewusr2)/(sigx4 *
      ewusr2) - sigx3 * ewusr2 * sigx9_2/(sigx4 * ewusr2)^2) *
      (1 - pwZ) * sigx2_2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1 * sigx2_1/ewusr1 -
      (sigx3 * sigx10/sigx4 + sigx1_2 * sigx2_2/ewusr2)) *
      dwZ/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr1) - 0.5) - 0.5 *
      (0.5 + 0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr1/ewusr1 -
      (sigx8_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv1/(2 *
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
    MARGIN = 1, STATS = wHvar * (((dmusig1 * ewvsr1/(4 *
      (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 *
      ewvsr1/ewusr1 + (0.5 * epsiu1 + 2 * sigx6_1) * sigx9_1))/(sigx4 *
      ewusr1) - (sigx8_1 * pwZ * sigx2_1 + 0.5 * (sigx4 *
      ewusr1)) * sigx9_1/(sigx4 * ewusr1)^2) * pwZ * sigx2_1,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * sigx8_1 * pwZ *
      ewusr2 * sigx2_1 * sigx2_2 * sigx9_2/((sigx4 * ewusr2)^2 *
      ewusr1)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx8_1 * (1 - sigx10 * pwZ/sigx4) * dwZ * sigx2_1/(sigx4 *
    ewusr1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr2 *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 + 0.5 * epsiu2 +
      2 * sigx6_2)) * dmusig2 * ewvsr2/ewusr2 - (sigx8_2 *
      (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 * epsiu2) *
      pmusig2))/(sigx4 * ewusr2) - ((1 - pwZ) * sigx8_2 *
      sigx2_2 + 0.5 * (sigx4 * ewusr2)) * sigx8_2/(sigx4 *
      ewusr2)^2) * (1 - pwZ) * sigx2_2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * sigx8_2 * pwZ *
      ewusr1 * sigx2_1 * sigx2_2 * sigx9_1/((sigx4 * ewusr1)^2 *
      ewusr2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((dmusig2 * ewvsr2/(4 * (ewu2 * ewusr2)) - 2 * (ewu2 *
      pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (epsiv2 *
      musig2) - 0.25) * dmusig2 * ewvsr2/ewusr2 + (0.5 *
      epsiu2 + 2 * sigx6_2) * sigx9_2))/(sigx4 * ewusr2) -
      ((1 - pwZ) * sigx8_2 * sigx2_2 + 0.5 * (sigx4 * ewusr2)) *
        sigx9_2/(sigx4 * ewusr2)^2) * (1 - pwZ) * sigx2_2,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (((1 - pwZ) * sigx10/sigx4 + 1) * sigx8_2 *
      dwZ * sigx2_2/(sigx4 * ewusr2)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewvsr1/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr1) - epsiv1^2 * musig1) *
      dmusig1)/(sigx4 * ewusr1) - pwZ * sigx2_1 * sigx9_1^2/(sigx4 *
      ewusr1)^2) * pwZ * sigx2_1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * ((1 - pwZ) * pwZ * ewusr2 * sigx2_1 *
      sigx2_2 * sigx9_1 * sigx9_2/((sigx4 * ewusr2)^2 *
      ewusr1)), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx10 * pwZ/sigx4) *
      dwZ * sigx2_1 * sigx9_1/(sigx4 * ewusr1), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_2/2 + (pmusig2 -
      epsiv2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewvsr2/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr2) - epsiv2^2 * musig2) *
      dmusig2)/(sigx4 * ewusr2) - (1 - pwZ) * sigx2_2 *
      sigx9_2^2/(sigx4 * ewusr2)^2) * (1 - pwZ) * sigx2_2,
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((1 - pwZ) * sigx10/sigx4 +
      1) * dwZ * sigx2_2 * sigx9_2/(sigx4 * ewusr2)), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar *
    (dwZ * (dwZ * sigx10/sigx4 + Wz) * sigx10/sigx4), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmcesfexponormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3 <- ((1 - prZ) * sigx1_1 * sigx2_1/ewusr1 + sigx1_2 *
    prZ * sigx2_2/ewusr2)
  sigx4 <- ((1 - prZ) * sigx2_1 * pmusig1/ewusr1 + prZ * sigx2_2 *
    pmusig2/ewusr2)
  sigx5_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx5_2 <- (dmusig2 * ewvsr2/ewusr2)
  epsiu1 <- (S * (epsilon)/ewusr1)
  epsiu2 <- (S * (epsilon)/ewusr2)
  sigx6_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx6_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx7_1 <- (0.5 * sigx5_1 - (0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx7_2 <- (0.5 * sigx5_2 - (0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  sigx8_1 <- (0.5 * sigx5_1 - (0.5 + 0.5 * epsiu1 + 2 * sigx6_1) *
    pmusig1)
  sigx8_2 <- (0.5 * sigx5_2 - (0.5 + 0.5 * epsiu2 + 2 * sigx6_2) *
    pmusig2)
  epsiv1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  epsiv2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx9_1 <- (ewv1 * pmusig1/(2 * ewu1) - epsiv1 * dmusig1)
  sigx9_2 <- (ewv2 * pmusig2/(2 * ewu2) - epsiv2 * dmusig2)
  sigx10 <- (sigx2_1 * pmusig1/ewusr1 - sigx2_2 * pmusig2/ewusr2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (((musig1/ewvsr1 - 1/ewusr1) * dmusig1/ewvsr1 -
      sigx1_1/ewusr1) * (1 - prZ) * sigx2_1/ewusr1 + ((musig2/ewvsr2 -
      1/ewusr2) * dmusig2/ewvsr2 - sigx1_2/ewusr2) * prZ *
      sigx2_2/ewusr2 - sigx3^2/sigx4)/sigx4, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu1 + 1 +
      2 * sigx6_1) * pmusig1 - 0.5 * sigx5_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 + 0.5 * epsiu1 + 2 *
        sigx6_1)/ewvsr1) * dmusig1)/(sigx4 * ewusr1) -
      sigx3 * sigx8_1 * ewusr1/(sigx4 * ewusr1)^2) * (1 -
      prZ) * sigx2_1, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * epsiu2 + 1 +
      2 * sigx6_2) * pmusig2 - 0.5 * sigx5_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 + 0.5 * epsiu2 + 2 *
        sigx6_2)/ewvsr2) * dmusig2)/(sigx4 * ewusr2) -
      sigx3 * sigx8_2 * ewusr2/(sigx4 * ewusr2)^2) * prZ *
      sigx2_2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((dmusig1 * (ewv1/(2 * ewu1) - (epsiv1 * musig1 +
    0.5))/ewvsr1 - sigx9_1/ewusr1)/(sigx4 * ewusr1) - sigx3 *
    ewusr1 * sigx9_1/(sigx4 * ewusr1)^2) * (1 - prZ) * sigx2_1,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((dmusig2 * (ewv2/(2 *
      ewu2) - (epsiv2 * musig2 + 0.5))/ewvsr2 - sigx9_2/ewusr2)/(sigx4 *
      ewusr2) - sigx3 * ewusr2 * sigx9_2/(sigx4 * ewusr2)^2) *
      prZ * sigx2_2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1 * sigx2_1/ewusr1 -
      (sigx3 * sigx10/sigx4 + sigx1_2 * sigx2_2/ewusr2)) *
      prZ * ewz/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr1) - 0.5) - 0.5 *
      (0.5 + 0.5 * epsiu1 + 2 * sigx6_1)) * dmusig1 * ewvsr1/ewusr1 -
      (sigx8_1 * (0.5 * epsiu1 + 2 * sigx6_1) + (2 * ((1 -
        8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv1/(2 *
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
    MARGIN = 1, STATS = wHvar * (((dmusig1 * ewvsr1/(4 *
      (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (epsiv1 * musig1) - 0.25) * dmusig1 *
      ewvsr1/ewusr1 + (0.5 * epsiu1 + 2 * sigx6_1) * sigx9_1))/(sigx4 *
      ewusr1) - (sigx8_1 * (1 - prZ) * sigx2_1 + 0.5 *
      (sigx4 * ewusr1)) * sigx9_1/(sigx4 * ewusr1)^2) *
      (1 - prZ) * sigx2_1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * sigx8_1 * (1 - prZ) *
      ewusr2 * sigx2_1 * sigx2_2 * sigx9_2/((sigx4 * ewusr2)^2 *
      ewusr1)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx8_1 * (1 - (1 - prZ) * sigx10/sigx4) * prZ * sigx2_1 *
    ewz/(sigx4 * ewusr1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (0.5 * (ewvsr2 *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 + 0.5 * epsiu2 +
      2 * sigx6_2)) * dmusig2 * ewvsr2/ewusr2 - (sigx8_2 *
      (0.5 * epsiu2 + 2 * sigx6_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 * epsiu2) *
      pmusig2))/(sigx4 * ewusr2) - (prZ * sigx8_2 * sigx2_2 +
      0.5 * (sigx4 * ewusr2)) * sigx8_2/(sigx4 * ewusr2)^2) *
      prZ * sigx2_2, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * sigx8_2 * (1 - prZ) *
      ewusr1 * sigx2_1 * sigx2_2 * sigx9_1/((sigx4 * ewusr1)^2 *
      ewusr2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((dmusig2 * ewvsr2/(4 * (ewu2 * ewusr2)) - 2 * (ewu2 *
      pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (epsiv2 *
      musig2) - 0.25) * dmusig2 * ewvsr2/ewusr2 + (0.5 *
      epsiu2 + 2 * sigx6_2) * sigx9_2))/(sigx4 * ewusr2) -
      (prZ * sigx8_2 * sigx2_2 + 0.5 * (sigx4 * ewusr2)) *
        sigx9_2/(sigx4 * ewusr2)^2) * prZ * sigx2_2,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx8_2 * (1 + prZ * sigx10/sigx4) *
      prZ * sigx2_2 * ewz/(sigx4 * ewusr2)), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_1/2 + (pmusig1 -
      epsiv1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewvsr1/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr1) - epsiv1^2 * musig1) *
      dmusig1)/(sigx4 * ewusr1) - (1 - prZ) * sigx2_1 *
      sigx9_1^2/(sigx4 * ewusr1)^2) * (1 - prZ) * sigx2_1,
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (prZ * (1 - prZ) * ewusr2 * sigx2_1 *
      sigx2_2 * sigx9_1 * sigx9_2/((sigx4 * ewusr2)^2 *
      ewusr1)), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 - prZ) * sigx10/sigx4) *
      prZ * sigx2_1 * sigx9_1 * ewz/(sigx4 * ewusr1), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx9_2/2 + (pmusig2 -
      epsiv2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewvsr2/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr2) - epsiv2^2 * musig2) *
      dmusig2)/(sigx4 * ewusr2) - prZ * sigx2_2 * sigx9_2^2/(sigx4 *
      ewusr2)^2) * prZ * sigx2_2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1 + prZ * sigx10/sigx4) *
      prZ * sigx2_2 * sigx9_2 * ewz/(sigx4 * ewusr2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    (1 - (1 + prZ * sigx10/sigx4) * ewz) * prZ * sigx10 *
    ewz/sigx4, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf exponential-normal distribution
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
# same sigma_u

## logit specification class membership
cnsfexponormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfexponormlike_logit(startVal, nXvar = nXvar,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfexponormlike_logit,
      grad = cgradcnsfexponormlike_logit, hess = chesscnsfexponormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfexponormlike_logit(mleObj$par,
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
      mleObj$hessian <- chesscnsfexponormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfexponormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfexponormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

## cauchit specification class membership
cnsfexponormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfexponormlike_cauchit(startVal, nXvar = nXvar,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfexponormlike_cauchit,
      grad = cgradcnsfexponormlike_cauchit, hess = chesscnsfexponormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfexponormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- chesscnsfexponormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfexponormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfexponormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

## probit specification class membership
cnsfexponormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfexponormlike_probit(startVal, nXvar = nXvar,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfexponormlike_probit,
      grad = cgradcnsfexponormlike_probit, hess = chesscnsfexponormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfexponormlike_probit(mleObj$par,
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
      mleObj$hessian <- chesscnsfexponormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfexponormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfexponormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

## cloglog specification class membership
cnsfexponormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfexponormlike_cloglog(startVal, nXvar = nXvar,
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
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfexponormlike_cloglog,
      grad = cgradcnsfexponormlike_cloglog, hess = chesscnsfexponormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfexponormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- chesscnsfexponormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfexponormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfexponormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

# different sigma_u

## logit specification class membership
mcesfexponormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfexponormlike_logit(startVal, nXvar = nXvar,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfexponormlike_logit,
      grad = cgradmcesfexponormlike_logit, hess = chessmcesfexponormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfexponormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfexponormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfexponormlike_logit(mleObj$par,
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
      mleObj$hessian <- chessmcesfexponormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfexponormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfexponormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfexponormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

## cauchit specification class membership
mcesfexponormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfexponormlike_cauchit(startVal, nXvar = nXvar,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfexponormlike_cauchit,
      grad = cgradmcesfexponormlike_cauchit, hess = chessmcesfexponormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfexponormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfexponormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfexponormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- chessmcesfexponormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfexponormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfexponormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfexponormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

## probit specification class membership
mcesfexponormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfexponormlike_probit(startVal, nXvar = nXvar,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfexponormlike_probit,
      grad = cgradmcesfexponormlike_probit, hess = chessmcesfexponormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfexponormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfexponormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfexponormlike_probit(mleObj$par,
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
      mleObj$hessian <- chessmcesfexponormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfexponormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfexponormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfexponormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

## cloglog specification class membership
mcesfexponormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initExpo <- start_st$initExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfexponormlike_cloglog(startVal, nXvar = nXvar,
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
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfexponormlike_cloglog,
      grad = cgradmcesfexponormlike_cloglog, hess = chessmcesfexponormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfexponormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfexponormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfexponormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- chessmcesfexponormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfexponormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfexponormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfexponormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initExpo = initExpo))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf exponential-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfexponormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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
}

## cauchit specification class membership
ccnsfexponormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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
}

## probit specification class membership
ccnsfexponormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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
}

## cloglog specification class membership
ccnsfexponormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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
}

# different sigma_u

## logit specification class membership
cmcesfexponormeff_logit <- function(object, level) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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

## cauchit specification class membership
cmcesfexponormeff_cauchit <- function(object, level) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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

## probit specification class membership
cmcesfexponormeff_probit <- function(object, level) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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

## cloglog specification class membership
cmcesfexponormeff_cloglog <- function(object, level) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sqrt(exp(Wv1)) * dnorm(mustar1/sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
  u_c2 <- mustar2 + sqrt(exp(Wv2)) * dnorm(mustar2/sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * exp(Wv1)) * pnorm(mustar1/sqrt(exp(Wv1)) -
      sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_c2 <- exp(-mustar2 + 1/2 * exp(Wv2)) * pnorm(mustar2/sqrt(exp(Wv2)) -
      sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * exp(Wv1)) *
      pnorm(mustar1/sqrt(exp(Wv1)) + sqrt(exp(Wv1)))/pnorm(mustar1/sqrt(exp(Wv1)))
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * exp(Wv2)) *
      pnorm(mustar2/sqrt(exp(Wv2)) + sqrt(exp(Wv2)))/pnorm(mustar2/sqrt(exp(Wv2)))
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
#' efficiencies for cnsf exponential-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfmargexponorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargexponorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
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
ccnsfmargexponorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargexponorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
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
ccnsfmargexponorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargexponorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
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
ccnsfmargexponorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargexponorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Pi1 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2))
  Pi2 <- 1/exp(Wu/2) * exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 *
    exp(Wu))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
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
cmcesfmargexponorm_Eu_logit <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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

cmcesfmargexponorm_Vu_logit <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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

## cauchit specification class membership
cmcesfmargexponorm_Eu_cauchit <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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

cmcesfmargexponorm_Vu_cauchit <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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

## probit specification class membership
cmcesfmargexponorm_Eu_probit <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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

cmcesfmargexponorm_Vu_probit <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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

## cloglog specification class membership
cmcesfmargexponorm_Eu_cloglog <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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

cmcesfmargexponorm_Vu_cloglog <- function(object) {
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
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
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
  mustar1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mustar2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Pi1 <- 1/exp(Wu1/2) * exp(object$S * epsilon/exp(Wu1/2) +
    exp(Wv1)/(2 * exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) -
    exp(Wv1/2)/exp(Wu1/2))
  Pi2 <- 1/exp(Wu2/2) * exp(object$S * epsilon/exp(Wu2/2) +
    exp(Wv2)/(2 * exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) -
    exp(Wv2/2)/exp(Wu2/2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
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
