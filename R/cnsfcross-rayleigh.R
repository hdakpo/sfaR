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
# Convolution: rayleigh - normal                                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf rayleigh-normal distribution
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
ccnsfraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
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
  mustar1 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# different sigma_u
cmcesfraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar1 <- -exp(Wu1) * S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf rayleigh-normal distribution
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
cstcnsfraynorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, itermax,
  printInfo, tol) {
  cat("Initialization: SFA + rayleigh - normal distributions...\n")
  initRay <- maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradraynormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initRay$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("CN_", colnames(Zvar)))
  names(initRay$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initRay = initRay))
}

# different sigma_u
cstmcesfraynorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + rayleigh - normal distributions...\n")
  initRay <- maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradraynormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = printInfo, reltol = tol), nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]),
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar)
  Esti <- initRay$estimate
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("MCE_",
    colnames(Zvar)))
  names(initRay$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initRay = initRay))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf rayleigh-normal distribution
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
cgradcnsfraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewu <- exp(Wu)
  ewz <- exp(Wz)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq1 <- ewu + ewv1
  sigma_sq2 <- ewu + ewv2
  sigmastar1 <- sqrt(ewu * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu * (epsilon)/ssq1)
  musig2 <- (S * ewu * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewvsq1 <- (1 - ewv1/(sigma_sq1))
  ewvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx1_1 <- (0.5 * (ewvsq1 * ewu/sigmastar1) + sigmastar1)
  sigx1_2 <- (0.5 * (ewvsq2 * ewu/sigmastar2) + sigmastar2)
  depsi1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv1))
  depsi2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv2))
  siguv1 <- (sigmastar1/ewv1 - ewu/ssq1)
  siguv2 <- (sigmastar2/ewv2 - ewu/ssq2)
  sigx2_1 <- (pmusig1 + S * dmusig1 * siguv1 * (epsilon))
  sigx2_2 <- (pmusig2 + S * dmusig2 * siguv2 * (epsilon))
  wusq1 <- (1 - ewu/(sigma_sq1))
  wusq2 <- (1 - ewu/(sigma_sq2))
  sigx3_1 <- (dmusig1 * sigmastar1 - S * ewu * pmusig1 * (epsilon)/(sigma_sq1))
  sigx3_2 <- (dmusig2 * sigmastar2 - S * ewu * pmusig2 * (epsilon)/(sigma_sq2))
  sigx4_1 <- (ewu * sigx2_1/(sigma_sq1) + S * wusq1 * sigx3_1 *
    (epsilon)/ewv1)
  sigx4_2 <- (ewu * sigx2_2/(sigma_sq2) + S * wusq2 * sigx3_2 *
    (epsilon)/ewv2)
  sigx5 <- (prC * depsi2 * sigx4_2 * sigmastar2/ewv2_h + depsi1 *
    sigx4_1 * ewz * sigmastar1/(wzdeno * ewv1_h))
  sigx6 <- (prC * sigx3_2 * depsi2 * sigmastar2/ewv2_h + sigx3_1 *
    depsi1 * ewz * sigmastar1/(wzdeno * ewv1_h))
  sigx7_1 <- (0.5 * (dmusig1 * ewv1/sigmastar1) - S * pmusig1 *
    (epsilon))
  sigx7_2 <- (0.5 * (dmusig2 * ewv2/sigmastar2) - S * pmusig2 *
    (epsilon))
  uvsig1 <- (wusq1 * ewv1/sigmastar1)
  uvsig2 <- (wusq2 * ewv2/sigmastar2)
  sigx8_1 <- (1/ssq1 - (0.5 * uvsig1 + sigmastar1) * ewu/ssq1^2)
  sigx8_2 <- (1/ssq2 - (0.5 * uvsig2 + sigmastar2) * ewu/ssq2^2)
  sigx9_1 <- (sigx7_1 * wusq1 + S^2 * sigx8_1 * sigx3_1 * ewu *
    (epsilon)^2/sigmastar1)
  sigx9_2 <- (sigx7_2 * wusq2 + S^2 * sigx8_2 * sigx3_2 * ewu *
    (epsilon)^2/sigmastar2)
  sigx10_1 <- (sigx9_1 * sigmastar1 + 0.5 * (wusq1 * sigx3_1 *
    ewv1/sigmastar1))
  sigx10_2 <- (sigx9_2 * sigmastar2 + 0.5 * (wusq2 * sigx3_2 *
    ewv2/sigmastar2))
  wzv1 <- (wzdeno * (sigma_sq1) * ewv1_h)
  wzv2 <- ((sigma_sq2) * ewv2_h)
  sigx11_1 <- (ewvsq1 * dmusig1/sigmastar1)
  sigx11_2 <- (ewvsq2 * dmusig2/sigmastar2)
  sigx12_1 <- (0.5 * sigx11_1 + S * pmusig1 * (epsilon)/(sigma_sq1))
  sigx12_2 <- (0.5 * sigx11_2 + S * pmusig2 * (epsilon)/(sigma_sq2))
  sigx13_1 <- (2 * ((S * (epsilon))^2/(2 * ewv1)^2) - S^2 *
    sigx1_1 * ewu^2 * (epsilon)^2/(ssq1^2 * (sigma_sq1) *
    sigmastar1))
  sigx13_2 <- (2 * ((S * (epsilon))^2/(2 * ewv2)^2) - S^2 *
    sigx1_2 * ewu^2 * (epsilon)^2/(ssq2^2 * (sigma_sq2) *
    sigmastar2))
  sigx14_1 <- (sigx12_1 * ewu/(sigma_sq1) + sigx13_1 * sigx3_1)
  sigx14_2 <- (sigx12_2 * ewu/(sigma_sq2) + sigx13_2 * sigx3_2)
  sigx15_1 <- (sigx14_1 * sigmastar1 + 0.5 * (ewvsq1 * sigx3_1 *
    ewu/ssq1))
  sigx15_2 <- (sigx14_2 * sigmastar2 + 0.5 * (ewvsq2 * sigx3_2 *
    ewu/ssq2))
  sigx16_1 <- (sigx15_1 * ewv1/(wzdeno * ewv1_h) - 0.5 * (wzdeno *
    sigx3_1 * ewv1_h * sigmastar1/(wzdeno * ewv1_h)^2))
  sigx16_2 <- (sigx15_2 * ewv2 - 0.5 * (sigx3_2 * sigmastar2))
  sigx17 <- (sigx10_1 * depsi1 * ewz/wzv1 + sigx10_2 * prC *
    depsi2/wzv2)
  sigx18 <- (1/(wzdeno * ewv1_h) - ewv1_h * ewz/(wzdeno * ewv1_h)^2)
  sigx19 <- (sigx18 * sigx3_1 * depsi1 * sigmastar1 - prC *
    sigx3_2 * depsi2 * sigmastar2/(wzdeno * ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx6,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx17 *
    ewu/sigx6 - 1), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = sigx16_1 * depsi1 * ewz/sigx6, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = sigx16_2 * prC * depsi2/(sigx6 *
      ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx19 *
    ewz/sigx6, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u
cgradmcesfraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu1 * (epsilon)/ssq1)
  musig2 <- (S * ewu2 * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv1))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv2))
  suvq1 <- (sigmastar1/ewv1 - ewu1/ssq1)
  suvq2 <- (sigmastar2/ewv2 - ewu2/ssq2)
  sigx2_1 <- (pmusig1 + S * dmusig1 * suvq1 * (epsilon))
  sigx2_2 <- (pmusig2 + S * dmusig2 * suvq2 * (epsilon))
  wusq1 <- (1 - ewu1/(sigma_sq1))
  wusq2 <- (1 - ewu2/(sigma_sq2))
  wvsq1 <- (1 - ewv1/(sigma_sq1))
  wvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx3_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/(sigma_sq1))
  sigx3_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/(sigma_sq2))
  uv1 <- (wzdeno * ewu1 * ewv1_h)
  uv2 <- (ewu2 * ewv2_h)
  sigx4_1 <- (ewu1 * sigx2_1/(sigma_sq1) + S * wusq1 * sigx3_1 *
    (epsilon)/ewv1)
  sigx4_2 <- (ewu2 * sigx2_2/(sigma_sq2) + S * wusq2 * sigx3_2 *
    (epsilon)/ewv2)
  sigx5 <- (prC * sigx1_2 * sigx4_2 * sigmastar2/uv2 + sigx1_1 *
    sigx4_1 * ewz * sigmastar1/uv1)
  sigx6 <- (prC * sigx3_2 * sigx1_2 * sigmastar2/uv2 + sigx3_1 *
    sigx1_1 * ewz * sigmastar1/uv1)
  dvsq1 <- (dmusig1 * ewv1/sigmastar1)
  dvsq2 <- (dmusig2 * ewv2/sigmastar2)
  sigx7_1 <- (0.5 * dvsq1 - S * pmusig1 * (epsilon))
  sigx7_2 <- (0.5 * dvsq2 - S * pmusig2 * (epsilon))
  wuvsq1 <- (wusq1 * ewv1/sigmastar1)
  wuvsq2 <- (wusq2 * ewv2/sigmastar2)
  s8sig1 <- (0.5 * wuvsq1 + sigmastar1)
  s8sig2 <- (0.5 * wuvsq2 + sigmastar2)
  sigx8_1 <- (1/ssq1 - s8sig1 * ewu1/ssq1^2)
  sigx8_2 <- (1/ssq2 - s8sig2 * ewu2/ssq2^2)
  sigx9_1 <- (sigx7_1 * wusq1 + S^2 * sigx8_1 * sigx3_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx9_2 <- (sigx7_2 * wusq2 + S^2 * sigx8_2 * sigx3_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx10_1 <- (wusq1 * sigx3_1 * ewv1/sigmastar1)
  sigx10_2 <- (wusq2 * sigx3_2 * ewv2/sigmastar2)
  sigx11_1 <- (sigx9_1 * sigmastar1 + 0.5 * sigx10_1)
  sigx11_2 <- (sigx9_2 * sigmastar2 + 0.5 * sigx10_2)
  sigx12_1 <- (wzdeno * (sigma_sq1) * ewv1_h)
  sigx12_2 <- ((sigma_sq2) * ewv2_h)
  sigx13_1 <- (sigx11_1/sigx12_1 - wzdeno * sigx3_1 * ewu1 *
    ewv1_h * sigmastar1/uv1^2)
  sigx13_2 <- (sigx11_2/sigx12_2 - sigx3_2 * ewu2 * ewv2_h *
    sigmastar2/uv2^2)
  sigx14_1 <- (wvsq1 * dmusig1/sigmastar1)
  sigx14_2 <- (wvsq2 * dmusig2/sigmastar2)
  sigx15_1 <- (0.5 * sigx14_1 + S * pmusig1 * (epsilon)/(sigma_sq1))
  sigx15_2 <- (0.5 * sigx14_2 + S * pmusig2 * (epsilon)/(sigma_sq2))
  sigx16_1 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx16_2 <- ((S * (epsilon))^2/(2 * ewv2)^2)
  sigx17_1 <- (wvsq1 * ewu1/sigmastar1)
  sigx17_2 <- (wvsq2 * ewu2/sigmastar2)
  sigx18_1 <- (0.5 * sigx17_1 + sigmastar1)
  sigx18_2 <- (0.5 * sigx17_2 + sigmastar2)
  sigx19_1 <- (ssq1^2 * (sigma_sq1) * sigmastar1)
  sigx19_2 <- (ssq2^2 * (sigma_sq2) * sigmastar2)
  sigx20_1 <- (2 * sigx16_1 - S^2 * sigx18_1 * ewu1^2 * (epsilon)^2/sigx19_1)
  sigx20_2 <- (2 * sigx16_2 - S^2 * sigx18_2 * ewu2^2 * (epsilon)^2/sigx19_2)
  sigx21_1 <- (sigx15_1 * ewu1/(sigma_sq1) + sigx20_1 * sigx3_1)
  sigx21_2 <- (sigx15_2 * ewu2/(sigma_sq2) + sigx20_2 * sigx3_2)
  sigx22_1 <- (wvsq1 * sigx3_1 * ewu1/ssq1)
  sigx22_2 <- (wvsq2 * sigx3_2 * ewu2/ssq2)
  sigx23_1 <- (wzdeno * sigx3_1 * ewu1 * ewv1_h * sigmastar1/uv1^2)
  sigx23_2 <- (sigx3_2 * ewu2 * ewv2_h * sigmastar2/uv2^2)
  sigx24_1 <- (sigx21_1 * sigmastar1 + 0.5 * sigx22_1)
  sigx24_2 <- (sigx21_2 * sigmastar2 + 0.5 * sigx22_2)
  sigx25_1 <- (sigx24_1 * ewv1/uv1 - 0.5 * sigx23_1)
  sigx25_2 <- (sigx24_2 * ewv2/uv2 - 0.5 * sigx23_2)
  sigx26_1 <- (1/uv1 - ewu1 * ewv1_h * ewz/uv1^2)
  sigx26_2 <- (wzdeno * ewu2 * ewv2_h)
  sigx27 <- (sigx26_1 * sigx3_1 * sigx1_1 * sigmastar1 - prC *
    sigx3_2 * sigx1_2 * sigmastar2/sigx26_2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx6,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13_1 *
    sigx1_1 * ewz/sigx6, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = sigx13_2 * prC * sigx1_2/sigx6, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = sigx25_1 * sigx1_1 * ewz/sigx6, FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx25_2 * prC * sigx1_2/sigx6,
      FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx27 *
      ewz/sigx6, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf rayleigh-normal distribution
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
chesscnsfraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewu <- exp(Wu)
  ewz <- exp(Wz)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  sigma_sq1 <- ewu + ewv1
  sigma_sq2 <- ewu + ewv2
  sigmastar1 <- sqrt(ewu * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu * (epsilon)/ssq1)
  musig2 <- (S * ewu * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewvsq1 <- (1 - ewv1/(sigma_sq1))
  ewvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx1_1 <- (0.5 * (ewvsq1 * ewu/sigmastar1) + sigmastar1)
  sigx1_2 <- (0.5 * (ewvsq2 * ewu/sigmastar2) + sigmastar2)
  depsi1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv1))
  depsi2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv2))
  siguv1 <- (sigmastar1/ewv1 - ewu/ssq1)
  siguv2 <- (sigmastar2/ewv2 - ewu/ssq2)
  sigx2_1 <- (pmusig1 + S * dmusig1 * siguv1 * (epsilon))
  sigx2_2 <- (pmusig2 + S * dmusig2 * siguv2 * (epsilon))
  wusq1 <- (1 - ewu/(sigma_sq1))
  wusq2 <- (1 - ewu/(sigma_sq2))
  sigx3_1 <- (dmusig1 * sigmastar1 - S * ewu * pmusig1 * (epsilon)/(sigma_sq1))
  sigx3_2 <- (dmusig2 * sigmastar2 - S * ewu * pmusig2 * (epsilon)/(sigma_sq2))
  sigx4_1 <- (ewu * sigx2_1/(sigma_sq1) + S * wusq1 * sigx3_1 *
    (epsilon)/ewv1)
  sigx4_2 <- (ewu * sigx2_2/(sigma_sq2) + S * wusq2 * sigx3_2 *
    (epsilon)/ewv2)
  sigx5 <- (prC * depsi2 * sigx4_2 * sigmastar2/ewv2_h + depsi1 *
    sigx4_1 * ewz * sigmastar1/(wzdeno * ewv1_h))
  sigx6 <- (prC * sigx3_2 * depsi2 * sigmastar2/ewv2_h + sigx3_1 *
    depsi1 * ewz * sigmastar1/(wzdeno * ewv1_h))
  sigx7_1 <- (0.5 * (dmusig1 * ewv1/sigmastar1) - S * pmusig1 *
    (epsilon))
  sigx7_2 <- (0.5 * (dmusig2 * ewv2/sigmastar2) - S * pmusig2 *
    (epsilon))
  uvsig1 <- (wusq1 * ewv1/sigmastar1)
  uvsig2 <- (wusq2 * ewv2/sigmastar2)
  sigx8_1 <- (1/ssq1 - (0.5 * uvsig1 + sigmastar1) * ewu/ssq1^2)
  sigx8_2 <- (1/ssq2 - (0.5 * uvsig2 + sigmastar2) * ewu/ssq2^2)
  sigx9_1 <- (sigx7_1 * wusq1 + S^2 * sigx8_1 * sigx3_1 * ewu *
    (epsilon)^2/sigmastar1)
  sigx9_2 <- (sigx7_2 * wusq2 + S^2 * sigx8_2 * sigx3_2 * ewu *
    (epsilon)^2/sigmastar2)
  sigx10_1 <- (sigx9_1 * sigmastar1 + 0.5 * (wusq1 * sigx3_1 *
    ewv1/sigmastar1))
  sigx10_2 <- (sigx9_2 * sigmastar2 + 0.5 * (wusq2 * sigx3_2 *
    ewv2/sigmastar2))
  wzv1 <- (wzdeno * (sigma_sq1) * ewv1_h)
  wzv2 <- ((sigma_sq2) * ewv2_h)
  sigx11_1 <- (ewvsq1 * dmusig1/sigmastar1)
  sigx11_2 <- (ewvsq2 * dmusig2/sigmastar2)
  sigx12_1 <- (0.5 * sigx11_1 + S * pmusig1 * (epsilon)/(sigma_sq1))
  sigx12_2 <- (0.5 * sigx11_2 + S * pmusig2 * (epsilon)/(sigma_sq2))
  sigx20_1 <- (ssq1^2 * (sigma_sq1) * sigmastar1)
  sigx20_2 <- (ssq2^2 * (sigma_sq2) * sigmastar2)
  sigx13_1 <- (2 * ((S * (epsilon))^2/(2 * ewv1)^2) - S^2 *
    sigx1_1 * ewu^2 * (epsilon)^2/sigx20_1)
  sigx13_2 <- (2 * ((S * (epsilon))^2/(2 * ewv2)^2) - S^2 *
    sigx1_2 * ewu^2 * (epsilon)^2/sigx20_2)
  sigx14_1 <- (sigx12_1 * ewu/(sigma_sq1) + sigx13_1 * sigx3_1)
  sigx14_2 <- (sigx12_2 * ewu/(sigma_sq2) + sigx13_2 * sigx3_2)
  sigx15_1 <- (sigx14_1 * sigmastar1 + 0.5 * (ewvsq1 * sigx3_1 *
    ewu/ssq1))
  sigx15_2 <- (sigx14_2 * sigmastar2 + 0.5 * (ewvsq2 * sigx3_2 *
    ewu/ssq2))
  sigx16_1 <- (sigx15_1 * ewv1/(wzdeno * ewv1_h) - 0.5 * (wzdeno *
    sigx3_1 * ewv1_h * sigmastar1/(wzdeno * ewv1_h)^2))
  sigx16_2 <- (sigx15_2 * ewv2 - 0.5 * (sigx3_2 * sigmastar2))
  sigx17 <- (sigx10_1 * depsi1 * ewz/wzv1 + sigx10_2 * prC *
    depsi2/wzv2)
  sigx18 <- (1/(wzdeno * ewv1_h) - ewv1_h * ewz/(wzdeno * ewv1_h)^2)
  sigx19 <- (sigx18 * sigx3_1 * depsi1 * sigmastar1 - prC *
    sigx3_2 * depsi2 * sigmastar2/(wzdeno * ewv2_h))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * ((((S^2 * ewu * (epsilon)^2/((sigma_sq1) *
      ewv1) - 1) * siguv1 + ewu/ssq1) * dmusig1 * ewu/(sigma_sq1) +
      wusq1 * (S * ((2 * (S * dmusig1 * siguv1 * (epsilon)) +
        3 * pmusig1) * ewu/(sigma_sq1) + S * wusq1 *
        sigx3_1 * (epsilon)/ewv1) * (epsilon) - dmusig1 *
        sigmastar1)/ewv1) * depsi1 * ewz * sigmastar1/(wzdeno *
      ewv1_h) + (((S^2 * ewu * (epsilon)^2/((sigma_sq2) *
      ewv2) - 1) * siguv2 + ewu/ssq2) * dmusig2 * ewu/(sigma_sq2) +
      wusq2 * (S * ((2 * (S * dmusig2 * siguv2 * (epsilon)) +
        3 * pmusig2) * ewu/(sigma_sq2) + S * wusq2 *
        sigx3_2 * (epsilon)/ewv2) * (epsilon) - dmusig2 *
        sigmastar2)/ewv2) * prC * depsi2 * sigmastar2/ewv2_h -
      sigx5^2/sigx6)/sigx6, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((wusq1 * (pmusig1 -
      0.5 * (S * dmusig1 * ewu * (epsilon)/ssq1)) + S *
      sigx8_1 * ewu * (S * ewu * sigx2_1 * (epsilon)/(sigma_sq1) -
      2 * sigx3_1) * (epsilon)/sigmastar1) * sigmastar1 +
      (0.5 * (ewu * ewv1 * sigx2_1/ssq1) + S * sigx10_1 *
        (epsilon)/ewv1) * wusq1) * depsi1 * ewz/wzv1 +
      ((wusq2 * (pmusig2 - 0.5 * (S * dmusig2 * ewu * (epsilon)/ssq2)) +
        S * sigx8_2 * ewu * (S * ewu * sigx2_2 * (epsilon)/(sigma_sq2) -
          2 * sigx3_2) * (epsilon)/sigmastar2) * sigmastar2 +
        (0.5 * (ewu * ewv2 * sigx2_2/ssq2) + S * sigx10_2 *
          (epsilon)/ewv2) * wusq2) * prC * depsi2/wzv2 -
      sigx17 * sigx5/sigx6) * ewu/sigx6, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (sigx16_1 * (S * wusq1 * (epsilon)/ewv1 - sigx5/sigx6) +
    (((sigx13_1 * sigx2_1 + (S * (0.5 * (ewvsq1/ewv1) + 1/(sigma_sq1)) *
      dmusig1 * ewu * (epsilon)/sigmastar1 - pmusig1)/(sigma_sq1)) *
      ewu/(sigma_sq1) + S * (2 * (sigx1_1 * ewu^2/sigx20_1) -
      4/(2 * ewv1)^2) * sigx3_1 * (epsilon)) * sigmastar1 +
      0.5 * (ewvsq1 * ewu^2 * sigx2_1/((sigma_sq1)^2 *
        sigmastar1))) * ewv1/(wzdeno * ewv1_h) - 0.5 *
    (wzdeno * ewu * ewv1_h * sigx2_1 * sigmastar1/((wzdeno *
      ewv1_h)^2 * (sigma_sq1)))) * depsi1 * ewz/sigx6,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((((sigx13_2 * sigx2_2 + (S * (0.5 *
      (ewvsq2/ewv2) + 1/(sigma_sq2)) * dmusig2 * ewu *
      (epsilon)/sigmastar2 - pmusig2)/(sigma_sq2)) * ewu/(sigma_sq2) +
      S * (2 * (sigx1_2 * ewu^2/sigx20_2) - 4/(2 * ewv2)^2) *
        sigx3_2 * (epsilon)) * sigmastar2 + 0.5 * (ewvsq2 *
      ewu^2 * sigx2_2/((sigma_sq2)^2 * sigmastar2))) *
      ewv2 + S * sigx16_2 * wusq2 * (epsilon)/ewv2 - 0.5 *
      (ewu * sigx2_2 * sigmastar2/(sigma_sq2)))/(sigx6 *
      ewv2_h) - sigx16_2 * sigx5 * ewv2_h/(sigx6 * ewv2_h)^2) *
      prC * depsi2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx18 * depsi1 * sigx4_1 *
      sigmastar1 - (sigx5 * sigx19/sigx6 + prC * depsi2 *
      sigx4_2 * sigmastar2/(wzdeno * ewv2_h))) * ewz/sigx6,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((((0.5 * ((sigx7_1 * wusq1 + S * ewu * pmusig1 * (epsilon)/(sigma_sq1) -
      dmusig1 * sigmastar1) * ewu/(sigma_sq1) - 0.5 * (wusq1 *
      sigx3_1)) + 0.5 * (sigx9_1 * ewu/(sigma_sq1))) *
      wusq1 * ewv1 + S^2 * sigx10_1 * sigx8_1 * ewu^2 *
      (epsilon)^2/(sigma_sq1))/sigmastar1 + (wusq1 * (ewu *
      (S^2 * sigx8_1 * dmusig1 * (epsilon)^2 - sigx7_1/(sigma_sq1)) -
      0.5 * (0.5 * (wusq1 * dmusig1 * ewv1/sigmastar1) +
        S^2 * sigx8_1 * dmusig1 * ewu * (epsilon)^2)) +
      S^2 * ((sigx7_1 * wusq1 * sigx8_1/(sigma_sq1) - ((0.5 *
        (ewu/(sigma_sq1)) + 1 - 0.5 * (0.5 * wusq1 +
        ewu/(sigma_sq1))) * wusq1 * ewv1/sigmastar1 +
        (2 - 2 * ((0.5 * uvsig1 + sigmastar1)^2 * ewu *
          (sigma_sq1)/ssq1^2)) * sigmastar1) * sigx3_1/ssq1^2) *
        ewu + (1 - 0.5 * wusq1) * sigx8_1 * sigx3_1) *
        ewu * (epsilon)^2/sigmastar1) * sigmastar1)/wzv1 +
      sigx10_1 * (1/wzv1 - wzdeno * ewu * ewv1_h/wzv1^2)) *
      depsi1 * ewz + ((((0.5 * ((sigx7_2 * wusq2 + S *
      ewu * pmusig2 * (epsilon)/(sigma_sq2) - dmusig2 *
      sigmastar2) * ewu/(sigma_sq2) - 0.5 * (wusq2 * sigx3_2)) +
      0.5 * (sigx9_2 * ewu/(sigma_sq2))) * wusq2 * ewv2 +
      S^2 * sigx10_2 * sigx8_2 * ewu^2 * (epsilon)^2/(sigma_sq2))/sigmastar2 +
      (wusq2 * (ewu * (S^2 * sigx8_2 * dmusig2 * (epsilon)^2 -
        sigx7_2/(sigma_sq2)) - 0.5 * (0.5 * (wusq2 *
        dmusig2 * ewv2/sigmastar2) + S^2 * sigx8_2 *
        dmusig2 * ewu * (epsilon)^2)) + S^2 * ((sigx7_2 *
        wusq2 * sigx8_2/(sigma_sq2) - ((0.5 * (ewu/(sigma_sq2)) +
        1 - 0.5 * (0.5 * wusq2 + ewu/(sigma_sq2))) *
        wusq2 * ewv2/sigmastar2 + (2 - 2 * ((0.5 * uvsig2 +
        sigmastar2)^2 * ewu * (sigma_sq2)/ssq2^2)) *
        sigmastar2) * sigx3_2/ssq2^2) * ewu + (1 - 0.5 *
        wusq2) * sigx8_2 * sigx3_2) * ewu * (epsilon)^2/sigmastar2) *
        sigmastar2)/wzv2 + sigx10_2 * (1/wzv2 - ewu *
      ewv2_h/wzv2^2)) * prC * depsi2 - sigx17^2 * ewu/sigx6) *
    ewu/sigx6, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (ewvsq1 * dmusig1) +
      0.5 * ((dmusig1 * ewv1/(sigma_sq1) - S^2 * ewvsq1 *
        sigx8_1 * dmusig1 * ewu * (epsilon)^2/sigmastar1) *
        ewu/(sigma_sq1) - 0.5 * (wusq1 * ewvsq1 * dmusig1)))/sigmastar1 +
      sigx7_1 * wusq1 * sigx13_1 + (S * (pmusig1 - ewu *
      (pmusig1/(sigma_sq1) + S * sigx8_1 * dmusig1 * (epsilon))) *
      (epsilon) - sigx12_1 * ewu)/(sigma_sq1))/(sigma_sq1) -
      S^2 * (((0.5 * (wusq1 * ewv1/(sigma_sq1)) + 0.5 *
        ((ewu/(sigma_sq1) - 1) * ewv1/(sigma_sq1) + 1 -
          0.5 * (wusq1 * ewvsq1))) * ewu/sigmastar1 +
        2 * sigx1_1)/sigx20_1 - ((ssq1^2 + 2 * ((0.5 *
        uvsig1 + sigmastar1) * (sigma_sq1)^2 * sigmastar1)) *
        sigmastar1 + 0.5 * (ssq1^2 * wusq1 * ewv1/sigmastar1)) *
        sigx1_1 * ewu/sigx20_1^2) * sigx3_1 * ewu * (epsilon)^2) *
      sigmastar1 + 0.5 * (((sigx7_1 * wusq1 * ewvsq1 +
      sigx3_1 * ewv1/(sigma_sq1)) * ewu/(sigma_sq1) + ewvsq1 *
      sigx3_1)/ssq1 - (0.5 * uvsig1 + sigmastar1) * ewvsq1 *
      sigx3_1 * ewu/ssq1^2) + 0.5 * (sigx14_1 * wusq1 *
      ewv1/ssq1)) * ewv1/(wzdeno * ewv1_h) + (S^2 * sigx16_1 *
      sigx8_1 * ewu * (epsilon)^2/sigmastar1 - 0.5 * ((sigx7_1 *
      sigmastar1 + 0.5 * (sigx3_1 * ewv1/sigmastar1)) *
      wusq1 * wzdeno * ewv1_h/(wzdeno * ewv1_h)^2))/(sigma_sq1) -
      sigx16_1 * sigx17/sigx6) * depsi1 * ewu * ewz/sigx6,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((0.5 * (ewvsq2 * dmusig2) +
      0.5 * ((dmusig2 * ewv2/(sigma_sq2) - S^2 * ewvsq2 *
        sigx8_2 * dmusig2 * ewu * (epsilon)^2/sigmastar2) *
        ewu/(sigma_sq2) - 0.5 * (wusq2 * ewvsq2 * dmusig2)))/sigmastar2 +
      sigx7_2 * wusq2 * sigx13_2 + (S * (pmusig2 - ewu *
      (pmusig2/(sigma_sq2) + S * sigx8_2 * dmusig2 * (epsilon))) *
      (epsilon) - sigx12_2 * ewu)/(sigma_sq2))/(sigma_sq2) -
      S^2 * (((0.5 * (wusq2 * ewv2/(sigma_sq2)) + 0.5 *
        ((ewu/(sigma_sq2) - 1) * ewv2/(sigma_sq2) + 1 -
          0.5 * (wusq2 * ewvsq2))) * ewu/sigmastar2 +
        2 * sigx1_2)/sigx20_2 - ((ssq2^2 + 2 * ((0.5 *
        uvsig2 + sigmastar2) * (sigma_sq2)^2 * sigmastar2)) *
        sigmastar2 + 0.5 * (ssq2^2 * wusq2 * ewv2/sigmastar2)) *
        sigx1_2 * ewu/sigx20_2^2) * sigx3_2 * ewu * (epsilon)^2) *
      sigmastar2 + 0.5 * (((sigx7_2 * wusq2 * ewvsq2 +
      sigx3_2 * ewv2/(sigma_sq2)) * ewu/(sigma_sq2) + ewvsq2 *
      sigx3_2)/ssq2 - (0.5 * uvsig2 + sigmastar2) * ewvsq2 *
      sigx3_2 * ewu/ssq2^2) + 0.5 * (sigx14_2 * wusq2 *
      ewv2/ssq2)) * ewv2 + (S^2 * sigx16_2 * sigx8_2 *
      ewu * (epsilon)^2/sigmastar2 - 0.5 * ((sigx7_2 *
      sigmastar2 + 0.5 * (sigx3_2 * ewv2/sigmastar2)) *
      wusq2))/(sigma_sq2))/(sigx6 * ewv2_h) - sigx16_2 *
      sigx17 * ewv2_h/(sigx6 * ewv2_h)^2) * prC * depsi2 *
      ewu, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx10_1 * sigx18 * depsi1/(sigma_sq1) -
      (sigx17 * sigx19/sigx6 + sigx10_2 * prC * depsi2/(wzdeno *
        (sigma_sq2) * ewv2_h))) * ewu * ewz/sigx6, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx12_1 * sigx13_1 +
      (S * (S * sigx1_1 * dmusig1 * ewu * (epsilon)/ssq1^2 -
        2 * (pmusig1/(sigma_sq1))) * (epsilon) - 0.5 *
        sigx11_1)/(sigma_sq1)) * ewv1 + (0.5 * (ewv1 *
      (S^2 * sigx1_1 * dmusig1 * ewu^2 * (epsilon)^2/(ssq1^2 *
        sigmastar1) - dmusig1)/(sigma_sq1) - 0.5 * (ewvsq1 *
      dmusig1)) + 0.5 * dmusig1) * ewvsq1/sigmastar1 +
      S * pmusig1 * (epsilon)/(sigma_sq1)) * ewu/(sigma_sq1) +
      ((2 - 16 * (ewv1^2/(2 * ewv1)^2)) * (S * (epsilon))^2/(2 *
        ewv1)^2 - S^2 * (sigx1_1 * (1/sigx20_1 - ((ssq1^2 +
        2 * (sigx1_1 * (sigma_sq1)^2 * sigmastar1)) *
        sigmastar1 + 0.5 * (ssq1^2 * ewvsq1 * ewu/sigmastar1)) *
        ewv1/sigx20_1^2) + (0.5 * (ewv1/(sigma_sq1)) -
        0.5 * (0.5 * ewvsq1 + ewv1/(sigma_sq1))) * ewvsq1/(ssq1^2 *
        ewv1)) * ewu^2 * (epsilon)^2) * sigx3_1) * sigmastar1 +
      ((0.5 * (((0.5 * sigx11_1 + 2 * (S * pmusig1 * (epsilon)/(sigma_sq1))) *
        ewu - dmusig1 * sigmastar1)/((sigma_sq1)^2 *
        sigmastar1) - sigx1_1 * sigx3_1/ssq1^2) + 0.5 *
        (sigx14_1/ssq1)) * ewv1 + 0.5 * (sigx3_1/ssq1)) *
        ewvsq1 * ewu)/(wzdeno * ewv1_h) + sigx16_1 *
      sigx13_1 - 0.5 * (sigx15_1 * wzdeno * ewv1_h/(wzdeno *
      ewv1_h)^2)) * ewv1 - (sigx16_1^2 * depsi1 * ewz/sigx6 +
      0.5 * (((sigx12_1 * ewu * ewv1/(sigma_sq1) + 0.5 *
        sigx3_1) * sigmastar1 + (0.5 * (ewvsq1 * ewu *
        ewv1/ssq1) - wzdeno^2 * ewv1_h^2 * sigmastar1/(wzdeno *
        ewv1_h)^2) * sigx3_1) * wzdeno * ewv1_h/(wzdeno *
        ewv1_h)^2))) * depsi1 * ewz/sigx6, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_1 * sigx16_2 * prC * depsi1 * depsi2 * ewv2_h *
      ewz/(sigx6 * ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx14_1 * sigx18 * ewv1 -
      ((0.5 - wzdeno^2 * ewv1_h^2/(wzdeno * ewv1_h)^2) *
        ewz + 0.5 * wzdeno) * sigx3_1 * ewv1_h/(wzdeno *
        ewv1_h)^2) * sigmastar1 + 0.5 * (ewvsq1 * sigx18 *
      sigx3_1 * ewu * ewv1/ssq1) - sigx16_1 * sigx19 *
      ewz/sigx6) * depsi1 * ewz/sigx6, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (((((sigx12_2 * sigx13_2 + (S * (S *
      sigx1_2 * dmusig2 * ewu * (epsilon)/ssq2^2 - 2 *
      (pmusig2/(sigma_sq2))) * (epsilon) - 0.5 * sigx11_2)/(sigma_sq2)) *
      ewv2 + (0.5 * (ewv2 * (S^2 * sigx1_2 * dmusig2 *
      ewu^2 * (epsilon)^2/(ssq2^2 * sigmastar2) - dmusig2)/(sigma_sq2) -
      0.5 * (ewvsq2 * dmusig2)) + 0.5 * dmusig2) * ewvsq2/sigmastar2 +
      S * pmusig2 * (epsilon)/(sigma_sq2)) * ewu/(sigma_sq2) +
      ((2 - 16 * (ewv2^2/(2 * ewv2)^2)) * (S * (epsilon))^2/(2 *
        ewv2)^2 - S^2 * (sigx1_2 * (1/sigx20_2 - ((ssq2^2 +
        2 * (sigx1_2 * (sigma_sq2)^2 * sigmastar2)) *
        sigmastar2 + 0.5 * (ssq2^2 * ewvsq2 * ewu/sigmastar2)) *
        ewv2/sigx20_2^2) + (0.5 * (ewv2/(sigma_sq2)) -
        0.5 * (0.5 * ewvsq2 + ewv2/(sigma_sq2))) * ewvsq2/(ssq2^2 *
        ewv2)) * ewu^2 * (epsilon)^2) * sigx3_2) * sigmastar2 +
      sigx16_2 * sigx13_2 + (((0.5 * (((0.5 * sigx11_2 +
      2 * (S * pmusig2 * (epsilon)/(sigma_sq2))) * ewu -
      dmusig2 * sigmastar2)/((sigma_sq2)^2 * sigmastar2) -
      sigx1_2 * sigx3_2/ssq2^2) + 0.5 * (sigx14_2/ssq2)) *
      ewv2 + 0.5 * (sigx3_2/ssq2)) * ewvsq2 - 0.5 * ((sigx12_2 *
      sigmastar2 + 0.5 * (ewvsq2 * sigx3_2/sigmastar2))/(sigma_sq2))) *
      ewu) * ewv2/(sigx6 * ewv2_h) - (sigx16_2 * prC *
      depsi2 + 0.5 * (sigx6 * ewv2_h)) * sigx16_2/(sigx6 *
      ewv2_h)^2) * prC * depsi2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx16_2 * sigx19/sigx6 +
      sigx15_2 * ewv2/wzdeno)/ewv2_h - 0.5 * (wzdeno *
      sigx3_2 * ewv2_h * sigmastar2/(wzdeno * ewv2_h)^2)) *
      prC * depsi2 * ewz/sigx6), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) +
      ewv2_h/(wzdeno * ewv2_h)^2) * sigx3_2 * depsi2 *
      sigmastar2 - (sigx19^2/sigx6 + (2 - 2 * (wzdeno *
      ewv1_h^2 * ewz/(wzdeno * ewv1_h)^2)) * sigx3_1 *
      depsi1 * ewv1_h * sigmastar1/(wzdeno * ewv1_h)^2)) *
      ewz + sigx18 * sigx3_1 * depsi1 * sigmastar1 - prC *
      sigx3_2 * depsi2 * sigmastar2/(wzdeno * ewv2_h)) *
      ewz/sigx6, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# different sigma_u
chessmcesfraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu1 * (epsilon)/ssq1)
  musig2 <- (S * ewu2 * (epsilon)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx1_1 <- exp(0.5 * (-musig1)^2 - (S * (epsilon))^2/(2 *
    ewv1))
  sigx1_2 <- exp(0.5 * (-musig2)^2 - (S * (epsilon))^2/(2 *
    ewv2))
  suvq1 <- (sigmastar1/ewv1 - ewu1/ssq1)
  suvq2 <- (sigmastar2/ewv2 - ewu2/ssq2)
  sigx2_1 <- (pmusig1 + S * dmusig1 * suvq1 * (epsilon))
  sigx2_2 <- (pmusig2 + S * dmusig2 * suvq2 * (epsilon))
  wusq1 <- (1 - ewu1/(sigma_sq1))
  wusq2 <- (1 - ewu2/(sigma_sq2))
  wvsq1 <- (1 - ewv1/(sigma_sq1))
  wvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx3_1 <- (dmusig1 * sigmastar1 - S * ewu1 * pmusig1 * (epsilon)/(sigma_sq1))
  sigx3_2 <- (dmusig2 * sigmastar2 - S * ewu2 * pmusig2 * (epsilon)/(sigma_sq2))
  uv1 <- (wzdeno * ewu1 * ewv1_h)
  uv2 <- (ewu2 * ewv2_h)
  sigx4_1 <- (ewu1 * sigx2_1/(sigma_sq1) + S * wusq1 * sigx3_1 *
    (epsilon)/ewv1)
  sigx4_2 <- (ewu2 * sigx2_2/(sigma_sq2) + S * wusq2 * sigx3_2 *
    (epsilon)/ewv2)
  sigx5 <- (prC * sigx1_2 * sigx4_2 * sigmastar2/uv2 + sigx1_1 *
    sigx4_1 * ewz * sigmastar1/uv1)
  sigx6 <- (prC * sigx3_2 * sigx1_2 * sigmastar2/uv2 + sigx3_1 *
    sigx1_1 * ewz * sigmastar1/uv1)
  dvsq1 <- (dmusig1 * ewv1/sigmastar1)
  dvsq2 <- (dmusig2 * ewv2/sigmastar2)
  sigx7_1 <- (0.5 * dvsq1 - S * pmusig1 * (epsilon))
  sigx7_2 <- (0.5 * dvsq2 - S * pmusig2 * (epsilon))
  wuvsq1 <- (wusq1 * ewv1/sigmastar1)
  wuvsq2 <- (wusq2 * ewv2/sigmastar2)
  s8sig1 <- (0.5 * wuvsq1 + sigmastar1)
  s8sig2 <- (0.5 * wuvsq2 + sigmastar2)
  sigx8_1 <- (1/ssq1 - s8sig1 * ewu1/ssq1^2)
  sigx8_2 <- (1/ssq2 - s8sig2 * ewu2/ssq2^2)
  sigx9_1 <- (sigx7_1 * wusq1 + S^2 * sigx8_1 * sigx3_1 * ewu1 *
    (epsilon)^2/sigmastar1)
  sigx9_2 <- (sigx7_2 * wusq2 + S^2 * sigx8_2 * sigx3_2 * ewu2 *
    (epsilon)^2/sigmastar2)
  sigx10_1 <- (wusq1 * sigx3_1 * ewv1/sigmastar1)
  sigx10_2 <- (wusq2 * sigx3_2 * ewv2/sigmastar2)
  sigx11_1 <- (sigx9_1 * sigmastar1 + 0.5 * sigx10_1)
  sigx11_2 <- (sigx9_2 * sigmastar2 + 0.5 * sigx10_2)
  sigx12_1 <- (wzdeno * (sigma_sq1) * ewv1_h)
  sigx12_2 <- ((sigma_sq2) * ewv2_h)
  sigx13_1 <- (sigx11_1/sigx12_1 - wzdeno * sigx3_1 * ewu1 *
    ewv1_h * sigmastar1/uv1^2)
  sigx13_2 <- (sigx11_2/sigx12_2 - sigx3_2 * ewu2 * ewv2_h *
    sigmastar2/uv2^2)
  sigx14_1 <- (wvsq1 * dmusig1/sigmastar1)
  sigx14_2 <- (wvsq2 * dmusig2/sigmastar2)
  sigx15_1 <- (0.5 * sigx14_1 + S * pmusig1 * (epsilon)/(sigma_sq1))
  sigx15_2 <- (0.5 * sigx14_2 + S * pmusig2 * (epsilon)/(sigma_sq2))
  sigx16_1 <- ((S * (epsilon))^2/(2 * ewv1)^2)
  sigx16_2 <- ((S * (epsilon))^2/(2 * ewv2)^2)
  sigx17_1 <- (wvsq1 * ewu1/sigmastar1)
  sigx17_2 <- (wvsq2 * ewu2/sigmastar2)
  sigx18_1 <- (0.5 * sigx17_1 + sigmastar1)
  sigx18_2 <- (0.5 * sigx17_2 + sigmastar2)
  sigx19_1 <- (ssq1^2 * (sigma_sq1) * sigmastar1)
  sigx19_2 <- (ssq2^2 * (sigma_sq2) * sigmastar2)
  sigx20_1 <- (2 * sigx16_1 - S^2 * sigx18_1 * ewu1^2 * (epsilon)^2/sigx19_1)
  sigx20_2 <- (2 * sigx16_2 - S^2 * sigx18_2 * ewu2^2 * (epsilon)^2/sigx19_2)
  sigx21_1 <- (sigx15_1 * ewu1/(sigma_sq1) + sigx20_1 * sigx3_1)
  sigx21_2 <- (sigx15_2 * ewu2/(sigma_sq2) + sigx20_2 * sigx3_2)
  sigx22_1 <- (wvsq1 * sigx3_1 * ewu1/ssq1)
  sigx22_2 <- (wvsq2 * sigx3_2 * ewu2/ssq2)
  sigx23_1 <- (wzdeno * sigx3_1 * ewu1 * ewv1_h * sigmastar1/uv1^2)
  sigx23_2 <- (sigx3_2 * ewu2 * ewv2_h * sigmastar2/uv2^2)
  sigx24_1 <- (sigx21_1 * sigmastar1 + 0.5 * sigx22_1)
  sigx24_2 <- (sigx21_2 * sigmastar2 + 0.5 * sigx22_2)
  sigx25_1 <- (sigx24_1 * ewv1/uv1 - 0.5 * sigx23_1)
  sigx25_2 <- (sigx24_2 * ewv2/uv2 - 0.5 * sigx23_2)
  sigx26_1 <- (1/uv1 - ewu1 * ewv1_h * ewz/uv1^2)
  sigx26_2 <- (wzdeno * ewu2 * ewv2_h)
  sigx27 <- (sigx26_1 * sigx3_1 * sigx1_1 * sigmastar1 - prC *
    sigx3_2 * sigx1_2 * sigmastar2/sigx26_2)
  sigx28_1 <- (S * wusq1 * (epsilon)/ewv1 - sigx5/sigx6)
  sigx28_2 <- (S * wusq2 * (epsilon)/ewv2 - sigx5/sigx6)
  sigx29_1 <- (sigx7_1 * wusq1 - S * pmusig1 * (epsilon))
  sigx29_2 <- (sigx7_2 * wusq2 - S * pmusig2 * (epsilon))
  sigx30_1 <- (sigx29_1 * ewu1/(sigma_sq1) + dmusig1 * sigmastar1)
  sigx30_2 <- (sigx29_2 * ewu2/(sigma_sq2) + dmusig2 * sigmastar2)
  sigx31_1 <- (sigx30_1 * sigmastar1 + (0.5 * (wusq1 * ewv1/ssq1) -
    2 * (wzdeno^2 * ewu1 * ewv1_h^2 * sigmastar1/uv1^2)) *
    sigx3_1 * ewu1)
  sigx31_2 <- (sigx30_2 * sigmastar2 + (0.5 * (wusq2 * ewv2/ssq2) -
    2 * (ewu2 * ewv2_h^2 * sigmastar2/uv2^2)) * sigx3_2 *
    ewu2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * ((((S^2 * ewu1 * (epsilon)^2/((sigma_sq1) *
      ewv1) - 1) * suvq1 + ewu1/ssq1) * dmusig1 * ewu1/(sigma_sq1) +
      wusq1 * (S * ((2 * (S * dmusig1 * suvq1 * (epsilon)) +
        3 * pmusig1) * ewu1/(sigma_sq1) + S * wusq1 *
        sigx3_1 * (epsilon)/ewv1) * (epsilon) - dmusig1 *
        sigmastar1)/ewv1) * sigx1_1 * ewz * sigmastar1/uv1 +
      (((S^2 * ewu2 * (epsilon)^2/((sigma_sq2) * ewv2) -
        1) * suvq2 + ewu2/ssq2) * dmusig2 * ewu2/(sigma_sq2) +
        wusq2 * (S * ((2 * (S * dmusig2 * suvq2 * (epsilon)) +
          3 * pmusig2) * ewu2/(sigma_sq2) + S * wusq2 *
          sigx3_2 * (epsilon)/ewv2) * (epsilon) - dmusig2 *
          sigmastar2)/ewv2) * prC * sigx1_2 * sigmastar2/uv2 -
      sigx5^2/sigx6)/sigx6, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx13_1 * sigx28_1 +
      (((wusq1 * (pmusig1 - 0.5 * (S * dmusig1 * ewu1 *
        (epsilon)/ssq1)) + S * sigx8_1 * ewu1 * (S *
        ewu1 * sigx2_1 * (epsilon)/(sigma_sq1) - 2 *
        sigx3_1) * (epsilon)/sigmastar1) * sigmastar1 +
        0.5 * (wusq1 * ewu1 * ewv1 * sigx2_1/ssq1))/(wzdeno *
        ewv1_h) - wzdeno * ewu1^2 * ewv1_h * sigx2_1 *
        sigmastar1/uv1^2)/(sigma_sq1)) * sigx1_1 * ewz/sigx6,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx13_2 * sigx28_2 +
      (((wusq2 * (pmusig2 - 0.5 * (S * dmusig2 * ewu2 *
        (epsilon)/ssq2)) + S * sigx8_2 * ewu2 * (S *
        ewu2 * sigx2_2 * (epsilon)/(sigma_sq2) - 2 *
        sigx3_2) * (epsilon)/sigmastar2) * sigmastar2 +
        0.5 * (wusq2 * ewu2 * ewv2 * sigx2_2/ssq2))/ewv2_h -
        ewu2^2 * ewv2_h * sigx2_2 * sigmastar2/uv2^2)/(sigma_sq2)) *
      prC * sigx1_2/sigx6, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (sigx25_1 * sigx28_1 + (((sigx20_1 * sigx2_1 + (S *
    (0.5 * (wvsq1/ewv1) + 1/(sigma_sq1)) * dmusig1 * ewu1 *
    (epsilon)/sigmastar1 - pmusig1)/(sigma_sq1)) * ewu1/(sigma_sq1) +
    S * (2 * (sigx18_1 * ewu1^2/sigx19_1) - 4/(2 * ewv1)^2) *
      sigx3_1 * (epsilon)) * sigmastar1 + 0.5 * (wvsq1 *
    ewu1^2 * sigx2_1/((sigma_sq1)^2 * sigmastar1))) * ewv1/uv1 -
    0.5 * (wzdeno * ewu1^2 * ewv1_h * sigx2_1 * sigmastar1/(uv1^2 *
      (sigma_sq1)))) * sigx1_1 * ewz/sigx6, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx25_2 * sigx28_2 +
      (((sigx20_2 * sigx2_2 + (S * (0.5 * (wvsq2/ewv2) +
        1/(sigma_sq2)) * dmusig2 * ewu2 * (epsilon)/sigmastar2 -
        pmusig2)/(sigma_sq2)) * ewu2/(sigma_sq2) + S *
        (2 * (sigx18_2 * ewu2^2/sigx19_2) - 4/(2 * ewv2)^2) *
        sigx3_2 * (epsilon)) * sigmastar2 + 0.5 * (wvsq2 *
        ewu2^2 * sigx2_2/((sigma_sq2)^2 * sigmastar2))) *
        ewv2/uv2 - 0.5 * (ewu2^2 * ewv2_h * sigx2_2 *
      sigmastar2/(uv2^2 * (sigma_sq2)))) * prC * sigx1_2/sigx6,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx26_1 * sigx1_1 *
      sigx4_1 * sigmastar1 - (sigx5 * sigx27/sigx6 + prC *
      sigx1_2 * sigx4_2 * sigmastar2/sigx26_2)) * ewz/sigx6,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((wusq1 * (ewu1 * (S^2 * sigx8_1 * dmusig1 * (epsilon)^2 -
      sigx7_1/(sigma_sq1)) - 0.5 * (0.5 * (wusq1 * dmusig1 *
      ewv1/sigmastar1) + S^2 * sigx8_1 * dmusig1 * ewu1 *
      (epsilon)^2)) + S^2 * ((sigx7_1 * wusq1 * sigx8_1/(sigma_sq1) -
      ((0.5 * (ewu1/(sigma_sq1)) + 1 - 0.5 * (0.5 * wusq1 +
        ewu1/(sigma_sq1))) * wusq1 * ewv1/sigmastar1 +
        (2 - 2 * (s8sig1^2 * ewu1 * (sigma_sq1)/ssq1^2)) *
          sigmastar1) * sigx3_1/ssq1^2) * ewu1 + (1 -
      0.5 * wusq1) * sigx8_1 * sigx3_1) * ewu1 * (epsilon)^2/sigmastar1) *
      sigmastar1 + (0.5 * ((sigx7_1 * wusq1 + S * ewu1 *
      pmusig1 * (epsilon)/(sigma_sq1) - dmusig1 * sigmastar1) *
      ewu1/(sigma_sq1) - 0.5 * (wusq1 * sigx3_1)) + 0.5 *
      (sigx9_1 * ewu1/(sigma_sq1))) * wusq1 * ewv1/sigmastar1)/sigx12_1 +
      ewu1 * (S^2 * sigx13_1 * sigx8_1 * ewu1 * (epsilon)^2/ssq1 -
        (sigx31_1/uv1^2 + sigx11_1/sigx12_1^2) * wzdeno *
          ewv1_h) - sigx13_1^2 * sigx1_1 * ewz/sigx6) *
    sigx1_1 * ewz/sigx6, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx13_1 * sigx13_2 * prC * sigx1_1 *
      sigx1_2 * ewz/sigx6^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((0.5 * (wvsq1 * dmusig1) +
      0.5 * ((dmusig1 * ewv1/(sigma_sq1) - S^2 * wvsq1 *
        sigx8_1 * dmusig1 * ewu1 * (epsilon)^2/sigmastar1) *
        ewu1/(sigma_sq1) - 0.5 * (wusq1 * wvsq1 * dmusig1)))/sigmastar1 +
      sigx7_1 * wusq1 * sigx20_1 + (S * (pmusig1 - ewu1 *
      (pmusig1/(sigma_sq1) + S * sigx8_1 * dmusig1 * (epsilon))) *
      (epsilon) - sigx15_1 * ewu1)/(sigma_sq1))/(sigma_sq1) -
      S^2 * (((0.5 * (wusq1 * ewv1/(sigma_sq1)) + 0.5 *
        ((ewu1/(sigma_sq1) - 1) * ewv1/(sigma_sq1) +
          1 - 0.5 * (wusq1 * wvsq1))) * ewu1/sigmastar1 +
        2 * sigx18_1)/sigx19_1 - ((ssq1^2 + 2 * (s8sig1 *
        (sigma_sq1)^2 * sigmastar1)) * sigmastar1 + 0.5 *
        (ssq1^2 * wusq1 * ewv1/sigmastar1)) * sigx18_1 *
        ewu1/sigx19_1^2) * sigx3_1 * ewu1 * (epsilon)^2) *
      sigmastar1 + 0.5 * (((sigx7_1 * wusq1 * wvsq1 + sigx3_1 *
      ewv1/(sigma_sq1)) * ewu1/(sigma_sq1) + wvsq1 * sigx3_1)/ssq1 -
      s8sig1 * wvsq1 * sigx3_1 * ewu1/ssq1^2) + 0.5 * (sigx21_1 *
      wusq1 * ewv1/ssq1))/(wzdeno * ewv1_h) - sigx24_1 *
      wzdeno * ewu1 * ewv1_h/uv1^2) * ewv1 + ewu1 * (S^2 *
      sigx25_1 * sigx8_1 * ewu1 * (epsilon)^2/ssq1 - 0.5 *
      (sigx31_1 * wzdeno * ewv1_h/uv1^2)) - sigx25_1 *
      sigx13_1 * sigx1_1 * ewz/sigx6) * sigx1_1 * ewz/sigx6,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx25_2 * sigx13_1 * prC *
      sigx1_1 * sigx1_2 * ewz/sigx6^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((sigx9_1 * sigx26_1/(sigma_sq1) - ((2 - 2 * (wzdeno^2 *
      ewu1^2 * ewv1_h^2/uv1^2)) * ewz + 1) * sigx3_1 *
      ewv1_h/uv1^2) * sigmastar1 + 0.5 * (wusq1 * sigx26_1 *
      sigx3_1 * ewv1/ssq1)) * ewu1 - sigx13_1 * sigx27 *
      ewz/sigx6) * sigx1_1 * ewz/sigx6, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((wusq2 * (ewu2 * (S^2 *
      sigx8_2 * dmusig2 * (epsilon)^2 - sigx7_2/(sigma_sq2)) -
      0.5 * (0.5 * (wusq2 * dmusig2 * ewv2/sigmastar2) +
        S^2 * sigx8_2 * dmusig2 * ewu2 * (epsilon)^2)) +
      S^2 * ((sigx7_2 * wusq2 * sigx8_2/(sigma_sq2) - ((0.5 *
        (ewu2/(sigma_sq2)) + 1 - 0.5 * (0.5 * wusq2 +
        ewu2/(sigma_sq2))) * wusq2 * ewv2/sigmastar2 +
        (2 - 2 * (s8sig2^2 * ewu2 * (sigma_sq2)/ssq2^2)) *
          sigmastar2) * sigx3_2/ssq2^2) * ewu2 + (1 -
        0.5 * wusq2) * sigx8_2 * sigx3_2) * ewu2 * (epsilon)^2/sigmastar2) *
      sigmastar2 + (0.5 * ((sigx7_2 * wusq2 + S * ewu2 *
      pmusig2 * (epsilon)/(sigma_sq2) - dmusig2 * sigmastar2) *
      ewu2/(sigma_sq2) - 0.5 * (wusq2 * sigx3_2)) + 0.5 *
      (sigx9_2 * ewu2/(sigma_sq2))) * wusq2 * ewv2/sigmastar2)/sigx12_2 +
      ewu2 * (S^2 * sigx13_2 * sigx8_2 * ewu2 * (epsilon)^2/ssq2 -
        (sigx31_2/uv2^2 + sigx11_2/sigx12_2^2) * ewv2_h) -
      sigx13_2^2 * prC * sigx1_2/sigx6) * prC * sigx1_2/sigx6,
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx25_1 * sigx13_2 * prC *
      sigx1_1 * sigx1_2 * ewz/sigx6^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((0.5 * (wvsq2 * dmusig2) + 0.5 * ((dmusig2 * ewv2/(sigma_sq2) -
      S^2 * wvsq2 * sigx8_2 * dmusig2 * ewu2 * (epsilon)^2/sigmastar2) *
      ewu2/(sigma_sq2) - 0.5 * (wusq2 * wvsq2 * dmusig2)))/sigmastar2 +
      sigx7_2 * wusq2 * sigx20_2 + (S * (pmusig2 - ewu2 *
      (pmusig2/(sigma_sq2) + S * sigx8_2 * dmusig2 * (epsilon))) *
      (epsilon) - sigx15_2 * ewu2)/(sigma_sq2))/(sigma_sq2) -
      S^2 * (((0.5 * (wusq2 * ewv2/(sigma_sq2)) + 0.5 *
        ((ewu2/(sigma_sq2) - 1) * ewv2/(sigma_sq2) +
          1 - 0.5 * (wusq2 * wvsq2))) * ewu2/sigmastar2 +
        2 * sigx18_2)/sigx19_2 - ((ssq2^2 + 2 * (s8sig2 *
        (sigma_sq2)^2 * sigmastar2)) * sigmastar2 + 0.5 *
        (ssq2^2 * wusq2 * ewv2/sigmastar2)) * sigx18_2 *
        ewu2/sigx19_2^2) * sigx3_2 * ewu2 * (epsilon)^2) *
      sigmastar2 + 0.5 * (((sigx7_2 * wusq2 * wvsq2 + sigx3_2 *
      ewv2/(sigma_sq2)) * ewu2/(sigma_sq2) + wvsq2 * sigx3_2)/ssq2 -
      s8sig2 * wvsq2 * sigx3_2 * ewu2/ssq2^2) + 0.5 * (sigx21_2 *
      wusq2 * ewv2/ssq2))/ewv2_h - sigx24_2 * ewu2 * ewv2_h/uv2^2) *
      ewv2 + ewu2 * (S^2 * sigx25_2 * sigx8_2 * ewu2 *
      (epsilon)^2/ssq2 - 0.5 * (sigx31_2 * ewv2_h/uv2^2)) -
      sigx25_2 * sigx13_2 * prC * sigx1_2/sigx6) * prC *
    sigx1_2/sigx6, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((sigx13_2 * sigx27/sigx6 + sigx11_2/(wzdeno *
      (sigma_sq2) * ewv2_h) - wzdeno * sigx3_2 * ewu2 *
      ewv2_h * sigmastar2/sigx26_2^2) * prC * sigx1_2 *
      ewz/sigx6), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx15_1 * sigx20_1 +
      (S * (S * sigx18_1 * dmusig1 * ewu1 * (epsilon)/ssq1^2 -
        2 * (pmusig1/(sigma_sq1))) * (epsilon) - 0.5 *
        sigx14_1)/(sigma_sq1)) * ewv1 + (0.5 * (ewv1 *
      (S^2 * sigx18_1 * dmusig1 * ewu1^2 * (epsilon)^2/(ssq1^2 *
        sigmastar1) - dmusig1)/(sigma_sq1) - 0.5 * (wvsq1 *
      dmusig1)) + 0.5 * dmusig1) * wvsq1/sigmastar1 + S *
      pmusig1 * (epsilon)/(sigma_sq1)) * ewu1/(sigma_sq1) +
      ((2 - 16 * (ewv1^2/(2 * ewv1)^2)) * (S * (epsilon))^2/(2 *
        ewv1)^2 - S^2 * (sigx18_1 * (1/sigx19_1 - ((ssq1^2 +
        2 * (sigx18_1 * (sigma_sq1)^2 * sigmastar1)) *
        sigmastar1 + 0.5 * (ssq1^2 * wvsq1 * ewu1/sigmastar1)) *
        ewv1/sigx19_1^2) + (0.5 * (ewv1/(sigma_sq1)) -
        0.5 * (0.5 * wvsq1 + ewv1/(sigma_sq1))) * wvsq1/(ssq1^2 *
        ewv1)) * ewu1^2 * (epsilon)^2) * sigx3_1) * sigmastar1 +
      ((0.5 * (((0.5 * sigx14_1 + 2 * (S * pmusig1 * (epsilon)/(sigma_sq1))) *
        ewu1 - dmusig1 * sigmastar1)/((sigma_sq1)^2 *
        sigmastar1) - sigx18_1 * sigx3_1/ssq1^2) + 0.5 *
        (sigx21_1/ssq1)) * ewv1 + 0.5 * (sigx3_1/ssq1)) *
        wvsq1 * ewu1)/uv1 + sigx25_1 * sigx20_1 - 0.5 *
      (sigx24_1 * wzdeno * ewu1 * ewv1_h/uv1^2)) * ewv1 -
      (sigx25_1^2 * sigx1_1 * ewz/sigx6 + 0.5 * (((sigx15_1 *
        ewu1 * ewv1/(sigma_sq1) + 0.5 * sigx3_1) * sigmastar1 +
        (0.5 * (wvsq1 * ewv1/ssq1) - wzdeno^2 * ewu1 *
          ewv1_h^2 * sigmastar1/uv1^2) * sigx3_1 * ewu1) *
        wzdeno * ewu1 * ewv1_h/uv1^2))) * sigx1_1 * ewz/sigx6,
    FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx25_1 * sigx25_2 * prC * sigx1_1 *
      sigx1_2 * ewz/sigx6^2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((sigx15_1 * sigx26_1 *
      ewv1/(sigma_sq1) - ((0.5 - wzdeno^2 * ewu1^2 * ewv1_h^2/uv1^2) *
      ewz + 0.5 * wzdeno) * sigx3_1 * ewv1_h/uv1^2) * ewu1 +
      sigx26_1 * sigx20_1 * sigx3_1 * ewv1) * sigmastar1 +
      0.5 * (wvsq1 * sigx26_1 * sigx3_1 * ewu1 * ewv1/ssq1) -
      sigx25_1 * sigx27 * ewz/sigx6) * sigx1_1 * ewz/sigx6,
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx15_2 * sigx20_2 +
      (S * (S * sigx18_2 * dmusig2 * ewu2 * (epsilon)/ssq2^2 -
        2 * (pmusig2/(sigma_sq2))) * (epsilon) - 0.5 *
        sigx14_2)/(sigma_sq2)) * ewv2 + (0.5 * (ewv2 *
      (S^2 * sigx18_2 * dmusig2 * ewu2^2 * (epsilon)^2/(ssq2^2 *
        sigmastar2) - dmusig2)/(sigma_sq2) - 0.5 * (wvsq2 *
      dmusig2)) + 0.5 * dmusig2) * wvsq2/sigmastar2 + S *
      pmusig2 * (epsilon)/(sigma_sq2)) * ewu2/(sigma_sq2) +
      ((2 - 16 * (ewv2^2/(2 * ewv2)^2)) * (S * (epsilon))^2/(2 *
        ewv2)^2 - S^2 * (sigx18_2 * (1/sigx19_2 - ((ssq2^2 +
        2 * (sigx18_2 * (sigma_sq2)^2 * sigmastar2)) *
        sigmastar2 + 0.5 * (ssq2^2 * wvsq2 * ewu2/sigmastar2)) *
        ewv2/sigx19_2^2) + (0.5 * (ewv2/(sigma_sq2)) -
        0.5 * (0.5 * wvsq2 + ewv2/(sigma_sq2))) * wvsq2/(ssq2^2 *
        ewv2)) * ewu2^2 * (epsilon)^2) * sigx3_2) * sigmastar2 +
      ((0.5 * (((0.5 * sigx14_2 + 2 * (S * pmusig2 * (epsilon)/(sigma_sq2))) *
        ewu2 - dmusig2 * sigmastar2)/((sigma_sq2)^2 *
        sigmastar2) - sigx18_2 * sigx3_2/ssq2^2) + 0.5 *
        (sigx21_2/ssq2)) * ewv2 + 0.5 * (sigx3_2/ssq2)) *
        wvsq2 * ewu2)/uv2 + sigx25_2 * sigx20_2 - 0.5 *
      (sigx24_2 * ewu2 * ewv2_h/uv2^2)) * ewv2 - (sigx25_2^2 *
      prC * sigx1_2/sigx6 + 0.5 * (((sigx15_2 * ewu2 *
      ewv2/(sigma_sq2) + 0.5 * sigx3_2) * sigmastar2 +
      (0.5 * (wvsq2 * ewv2/ssq2) - ewu2 * ewv2_h^2 * sigmastar2/uv2^2) *
        sigx3_2 * ewu2) * ewu2 * ewv2_h/uv2^2))) * prC *
      sigx1_2/sigx6, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx25_2 * sigx27/sigx6 +
      sigx24_2 * ewv2/sigx26_2 - 0.5 * (wzdeno * sigx3_2 *
      ewu2 * ewv2_h * sigmastar2/sigx26_2^2)) * prC * sigx1_2 *
      ewz/sigx6), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    ((prC * (1/(wzdeno^2 * ewu2 * ewv2_h) + ewu2 * ewv2_h/sigx26_2^2) *
      sigx3_2 * sigx1_2 * sigmastar2 - (sigx27^2/sigx6 +
      (2 - 2 * (wzdeno * ewu1^2 * ewv1_h^2 * ewz/uv1^2)) *
        sigx3_1 * sigx1_1 * ewu1 * ewv1_h * sigmastar1/uv1^2)) *
      ewz + sigx26_1 * sigx3_1 * sigx1_1 * sigmastar1 -
      prC * sigx3_2 * sigx1_2 * sigmastar2/sigx26_2) *
    ewz/sigx6, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf rayleigh-normal distribution
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
cnsfraynormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfraynormlike(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(ccnsfraynormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfraynormlike, grad = cgradcnsfraynormlike,
      hess = chesscnsfraynormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfraynormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfraynormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfraynormlike(mleObj$par,
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
      mleObj$hessian <- chesscnsfraynormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfraynormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfraynormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfraynormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initRay = initRay))
}

# different sigma_u
mcesfraynormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initRay <- start_st$initRay
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfraynormlike(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmcesfraynormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfraynormlike, grad = cgradmcesfraynormlike,
      hess = chessmcesfraynormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfraynormlike(mleObj$par,
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
      mleObj$hessian <- chessmcesfraynormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfraynormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfraynormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfraynormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initRay = initRay))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf rayleigh-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u
ccnsfraynormeff <- function(object, level) {
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu)) * sigmastar2 * (sigmastar2 *
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

# different sigma_u
cmcesfraynormeff <- function(object, level) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
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
#' marginal impact on efficiencies for cnsf rayleigh-normal distribution
#' @param object object of class sfacross
#' @noRd
# same sigma_u
ccnsfmargraynorm_Eu <- function(object) {
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * 1/2 * sqrt(pi/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmargraynorm_Vu <- function(object) {
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
  mustar1 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu)) * sigmastar2 * (sigmastar2 *
    dnorm(mustar2/sigmastar2) + mustar2 * pnorm(mustar2/sigmastar2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (4 - pi)/2, ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (4 - pi)/2, ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# different sigma_u
cmcesfmargraynorm_Eu <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
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
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmargraynorm_Vu <- function(object) {
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
  mustar1 <- -exp(Wu1) * object$S * epsilon/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- exp(1/2 * (mustar1/sigmastar1)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv1)))/(exp(Wv1/2) * exp(Wu1)) * sigmastar1 * (sigmastar1 *
    dnorm(mustar1/sigmastar1) + mustar1 * pnorm(mustar1/sigmastar1))
  Pi2 <- exp(1/2 * (mustar2/sigmastar2)^2 - (object$S * epsilon)^2/(2 *
    exp(Wv2)))/(exp(Wv2/2) * exp(Wu2)) * sigmastar2 * (sigmastar2 *
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
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
