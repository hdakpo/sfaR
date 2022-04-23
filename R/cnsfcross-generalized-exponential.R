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
# Convolution: generalized exponential - normal                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf generalized exponential-normal distribution
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
ccnsfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  A1 <- S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B1 <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a1 <- -S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b1 <- -S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  A2 <- S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))
  B2 <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv2)/exp(Wu)
  a2 <- -S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)
  b2 <- -S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# different sigma_u
cmcesfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  A1 <- S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))
  B1 <- 2 * S * epsilon/exp(Wu1/2) + 2 * exp(Wv1)/exp(Wu1)
  a1 <- -S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)
  b1 <- -S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu1/2)
  A2 <- S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))
  B2 <- 2 * S * epsilon/exp(Wu2/2) + 2 * exp(Wv2)/exp(Wu2)
  a2 <- -S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)
  b2 <- -S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf generalized exponential-normal distribution
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
cstcnsfgenexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + generalized exponential - normal distributions...\n")
  initGenExpo <- maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradgenexponormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initGenExpo$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("CN_", colnames(Zvar)))
  names(initGenExpo$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# different sigma_u
cstmcesfgenexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + generalized exponential - normal distributions...\n")
  initGenExpo <- maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradgenexponormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = printInfo,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initGenExpo$estimate
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("MCE_",
    colnames(Zvar)))
  names(initGenExpo$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf generalized exponential-normal distribution
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
cgradcnsfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  musig1 <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewv1_h/ewu_h) + S * (epsilon)/ewv1_h)
  epsivu2 <- (2 * (ewv2_h/ewu_h) + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu_h)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx9 <- (S * (epsilon)/ewu_h)
  sigx3_1 <- exp(2 * (ewv1/ewu) + 2 * sigx9)
  sigx3_2 <- exp(2 * (ewv2/ewu) + 2 * sigx9)
  sigx4_1 <- (sigx1_1 * sigx2_1 - (depsi1/ewv1_h - 2 * (pepsi1/ewu_h)) *
    sigx3_1)
  sigx4_2 <- (sigx1_2 * sigx2_2 - (depsi2/ewv2_h - 2 * (pepsi2/ewu_h)) *
    sigx3_2)
  sigx5 <- (2 * (sigx4_1 * ewz/wzdeno) + 2 * (sigx4_2 * prC))
  sigx6_1 <- (sigx2_1 * pmusig1 - sigx3_1 * pepsi1)
  sigx19 <- (sigx2_2 * pmusig2 - sigx3_2 * pepsi2)
  sigx6_2 <- (prC * sigx19)
  sigx7 <- (2 * sigx6_2 + 2 * (sigx6_1 * ewz/wzdeno))
  sigx8_1 <- (dmusig1 * ewv1_h/ewu_h)
  sigx8_2 <- (dmusig2 * ewv2_h/ewu_h)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- (2 * (ewv1/ewu) + S * (epsilon)/ewu_h)
  sigx12_2 <- (2 * (ewv2/ewu) + S * (epsilon)/ewu_h)
  sigx13_1 <- (depsi1 * ewv1_h/ewu_h - sigx12_1 * pepsi1)
  sigx13_2 <- (depsi2 * ewv2_h/ewu_h - sigx12_2 * pepsi2)
  sigx14_1 <- ((sigx11_1 * sigx2_1 - sigx13_1 * sigx3_1) *
    ewz/wzdeno)
  sigx14_2 <- ((sigx11_2 * sigx2_2 - sigx13_2 * sigx3_2) *
    prC)
  sigx15_1 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx15_2 <- (0.5 * (ewv2_h/ewu_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (ewv1_h/ewu_h - 0.5 * (S * (epsilon)/ewv1_h))
  sigx17_2 <- (ewv2_h/ewu_h - 0.5 * (S * (epsilon)/ewv2_h))
  sigx18_1 <- (2 * (ewv1 * pepsi1/ewu) - depsi1 * sigx17_1)
  sigx18_2 <- (2 * (ewv2 * pepsi2/ewu) - depsi2 * sigx17_2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ((2 * sigx14_1 +
    2 * sigx14_2)/sigx7 - 0.5), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * ((sigx2_1 * sigx16_1 - sigx18_1 *
      sigx3_1) * ewz/(wzdeno * sigx7)), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (prC * (sigx2_2 * sigx16_2 -
      sigx18_2 * sigx3_2)/sigx7), FUN = "*"), sweep(Zvar,
    MARGIN = 1, STATS = prC * (2 * sigx6_1 - 2 * sigx19) *
      ewz/(wzdeno * sigx7), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u
cgradmcesfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  epsi1 <- (2 * (ewv1_h/ewu1_h) + S * (epsilon)/ewv1_h)
  epsi2 <- (2 * (ewv2_h/ewu2_h) + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsi1, 0, 1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  pepsi2 <- pnorm(-epsi2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu1_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu2_h)
  sigx2_1 <- (depsi1/ewv1_h - 2 * (pepsi1/ewu1_h))
  sigx2_2 <- (depsi2/ewv2_h - 2 * (pepsi2/ewu2_h))
  sigx3_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx3_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx4_1 <- exp(2 * (ewv1/ewu1) + 2 * (S * (epsilon)/ewu1_h))
  sigx4_2 <- exp(2 * (ewv2/ewu2) + 2 * (S * (epsilon)/ewu2_h))
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
  sigx10_1 <- (dmusig1 * ewv1_h/ewu1_h)
  sigx10_2 <- (dmusig2 * ewv2_h/ewu2_h)
  sigx11_1 <- (0.5 * (S * (epsilon)/ewu1_h) + 2 * (ewu1 * ewv1/(2 *
    ewu1)^2))
  sigx11_2 <- (0.5 * (S * (epsilon)/ewu2_h) + 2 * (ewu2 * ewv2/(2 *
    ewu2)^2))
  uvsi1 <- (2 * (ewv1/ewu1) + S * (epsilon)/ewu1_h)
  uvsi2 <- (2 * (ewv2/ewu2) + S * (epsilon)/ewu2_h)
  sigx12_1 <- (depsi1 * ewv1_h/ewu1_h - uvsi1 * pepsi1)
  sigx12_2 <- (depsi2 * ewv2_h/ewu2_h - uvsi2 * pepsi2)
  sp1 <- (0.5 * sigx10_1 - sigx11_1 * pmusig1)
  sp2 <- (0.5 * sigx10_2 - sigx11_2 * pmusig2)
  sigx13_1 <- (sp1 * sigx3_1 - sigx12_1 * sigx4_1)
  sigx13_2 <- (sp2 * sigx3_2 - (sigx12_2 * sigx4_2 + 0.5 *
    sigx7_2))
  sigx14_1 <- (wzdeno * ewu1_h * sigx7_1/(wzdeno * ewu1_h)^2)
  sigx15_1 <- (sigx13_1/(wzdeno * ewu1_h) - 0.5 * sigx14_1)
  sigx16_1 <- (0.5 * (ewv1_h/ewu1_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx16_2 <- (0.5 * (ewv2_h/ewu2_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx17_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx16_1 * dmusig1)
  sigx17_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx16_2 * dmusig2)
  sigx18_1 <- (ewv1_h/ewu1_h - 0.5 * (S * (epsilon)/ewv1_h))
  sigx18_2 <- (ewv2_h/ewu2_h - 0.5 * (S * (epsilon)/ewv2_h))
  sigx19_1 <- (sigx3_1 * sigx17_1 - (2 * (ewv1 * pepsi1/ewu1) -
    depsi1 * sigx18_1) * sigx4_1)
  sigx19_2 <- (sigx3_2 * sigx17_2 - (2 * (ewv2 * pepsi2/ewu2) -
    depsi2 * sigx18_2) * sigx4_2)
  sigx20_1 <- (1/(wzdeno * ewu1_h) - ewu1_h * ewz/(wzdeno *
    ewu1_h)^2)
  sigx20_2 <- (prC * sigx7_2/(wzdeno * ewu2_h))
  sigx21 <- (2 * (sigx20_1 * sigx7_1) - 2 * sigx20_2)
  wzsig <- (wzdeno * s8 * ewu1_h)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx9,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx15_1 *
    ewz/s8), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx13_2 * prC/(s8 * ewu2_h)), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (sigx19_1 * ewz/wzsig), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (prC * sigx19_2/(s8 *
      ewu2_h)), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx21 *
      ewz/s8, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf generalized exponential-normal distribution
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
chesscnsfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  musig1 <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- (2 * (ewv1_h/ewu_h) + S * (epsilon)/ewv1_h)
  epsivu2 <- (2 * (ewv2_h/ewu_h) + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu_h)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx9 <- (S * (epsilon)/ewu_h)
  sigx3_1 <- exp(2 * (ewv1/ewu) + 2 * sigx9)
  sigx3_2 <- exp(2 * (ewv2/ewu) + 2 * sigx9)
  sigx4_1 <- (sigx1_1 * sigx2_1 - (depsi1/ewv1_h - 2 * (pepsi1/ewu_h)) *
    sigx3_1)
  sigx4_2 <- (sigx1_2 * sigx2_2 - (depsi2/ewv2_h - 2 * (pepsi2/ewu_h)) *
    sigx3_2)
  sigx5 <- (2 * (sigx4_1 * ewz/wzdeno) + 2 * (sigx4_2 * prC))
  sigx6_1 <- (sigx2_1 * pmusig1 - sigx3_1 * pepsi1)
  sigx19 <- (sigx2_2 * pmusig2 - sigx3_2 * pepsi2)
  sigx6_2 <- (prC * sigx19)
  sigx7 <- (2 * sigx6_2 + 2 * (sigx6_1 * ewz/wzdeno))
  sigx8_1 <- (dmusig1 * ewv1_h/ewu_h)
  sigx8_2 <- (dmusig2 * ewv2_h/ewu_h)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- (2 * (ewv1/ewu) + S * (epsilon)/ewu_h)
  sigx12_2 <- (2 * (ewv2/ewu) + S * (epsilon)/ewu_h)
  sigx13_1 <- (depsi1 * ewv1_h/ewu_h - sigx12_1 * pepsi1)
  sigx13_2 <- (depsi2 * ewv2_h/ewu_h - sigx12_2 * pepsi2)
  sigx14_1 <- ((sigx11_1 * sigx2_1 - sigx13_1 * sigx3_1) *
    ewz/wzdeno)
  sigx14_2 <- ((sigx11_2 * sigx2_2 - sigx13_2 * sigx3_2) *
    prC)
  sigx15_1 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx15_2 <- (0.5 * (ewv2_h/ewu_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (ewv1_h/ewu_h - 0.5 * (S * (epsilon)/ewv1_h))
  sigx17_2 <- (ewv2_h/ewu_h - 0.5 * (S * (epsilon)/ewv2_h))
  sigx18_1 <- (2 * (ewv1 * pepsi1/ewu) - depsi1 * sigx17_1)
  sigx18_2 <- (2 * (ewv2 * pepsi2/ewu) - depsi2 * sigx17_2)
  ewusq <- (1 - 8 * (ewu^2/(2 * ewu)^2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (2 * ((((musig1/ewv1_h - 1/ewu_h) *
      dmusig1/ewv1_h - sigx1_1/ewu_h) * sigx2_1 - ((epsivu1/ewv1_h -
      2/ewu_h) * depsi1/ewv1_h - 2 * ((depsi1/ewv1_h -
      2 * (pepsi1/ewu_h))/ewu_h)) * sigx3_1) * ewz/wzdeno) +
      2 * ((((musig2/ewv2_h - 1/ewu_h) * dmusig2/ewv2_h -
        sigx1_2/ewu_h) * sigx2_2 - ((epsivu2/ewv2_h -
        2/ewu_h) * depsi2/ewv2_h - 2 * ((depsi2/ewv2_h -
        2 * (pepsi2/ewu_h))/ewu_h)) * sigx3_2) * prC) -
      sigx5^2/sigx7)/sigx7, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((((0.5 + 0.5 *
      sigx9 + 2 * sigx10_1) * pmusig1 - 0.5 * sigx8_1)/ewu_h +
      (0.5 * (musig1/ewu_h) - (0.5 * sigx9 + 2 * sigx10_1)/ewv1_h) *
        dmusig1) * sigx2_1 - ((epsivu1/ewu_h - sigx12_1/ewv1_h) *
      depsi1 + (pepsi1 - 2 * sigx13_1)/ewu_h) * sigx3_1) *
      ewz/wzdeno) + 2 * (((((0.5 + 0.5 * sigx9 + 2 * sigx10_2) *
      pmusig2 - 0.5 * sigx8_2)/ewu_h + (0.5 * (musig2/ewu_h) -
      (0.5 * sigx9 + 2 * sigx10_2)/ewv2_h) * dmusig2) *
      sigx2_2 - ((epsivu2/ewu_h - sigx12_2/ewv2_h) * depsi2 +
      (pepsi2 - 2 * sigx13_2)/ewu_h) * sigx3_2) * prC) -
      (2 * sigx14_1 + 2 * sigx14_2) * sigx5/sigx7)/sigx7,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((dmusig1 * (ewv1/(2 * ewu) - (sigx15_1 * musig1 +
    0.5))/ewv1_h - sigx16_1/ewu_h) * sigx2_1 - ((2 * (ewv1/ewu) -
    (epsivu1 * sigx17_1 + 0.5)) * depsi1/ewv1_h - 2 * (sigx18_1/ewu_h)) *
    sigx3_1)/(wzdeno * sigx7) - wzdeno * sigx5 * (sigx2_1 *
    sigx16_1 - sigx18_1 * sigx3_1)/(wzdeno * sigx7)^2) *
    ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * ((dmusig2 * (ewv2/(2 * ewu) -
      (sigx15_2 * musig2 + 0.5))/ewv2_h - sigx16_2/ewu_h) *
      sigx2_2 - (((2 * (ewv2/ewu) - (epsivu2 * sigx17_2 +
      0.5)) * depsi2/ewv2_h - 2 * (sigx18_2/ewu_h)) * sigx3_2 +
      sigx5 * (sigx2_2 * sigx16_2 - sigx18_2 * sigx3_2)/sigx7)) *
      prC/sigx7), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * sigx4_1 - 2 * sigx4_2)/(wzdeno *
      sigx7) - wzdeno * sigx5 * (2 * sigx6_1 - 2 * sigx19)/(wzdeno *
      sigx7)^2) * prC * ewz, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((((0.5 * (0.5 * (ewv1_h * musig1/ewu_h) - 0.5) -
      0.5 * (0.5 * sigx9 + 2 * sigx10_1)) * dmusig1 * ewv1_h/ewu_h -
      (sigx11_1 * (0.5 * sigx9 + 2 * sigx10_1) + (2 * (ewusq *
        ewu * ewv1/(2 * ewu)^2) - 0.25 * sigx9) * pmusig1)) *
      sigx2_1 - (((epsivu1 * ewv1_h - S * (epsilon))/ewu_h -
      (0.5 + 2 * (ewv1/ewu))) * depsi1 * ewv1_h/ewu_h +
      (0.5 * sigx9 + 2 * (ewv1/ewu)) * pepsi1 - sigx12_1 *
      sigx13_1) * sigx3_1) * ewz/wzdeno) + 2 * ((((0.5 *
      (0.5 * (ewv2_h * musig2/ewu_h) - 0.5) - 0.5 * (0.5 *
      sigx9 + 2 * sigx10_2)) * dmusig2 * ewv2_h/ewu_h -
      (sigx11_2 * (0.5 * sigx9 + 2 * sigx10_2) + (2 * (ewusq *
        ewu * ewv2/(2 * ewu)^2) - 0.25 * sigx9) * pmusig2)) *
      sigx2_2 - (((epsivu2 * ewv2_h - S * (epsilon))/ewu_h -
      (0.5 + 2 * (ewv2/ewu))) * depsi2 * ewv2_h/ewu_h +
      (0.5 * sigx9 + 2 * (ewv2/ewu)) * pepsi2 - sigx12_2 *
      sigx13_2) * sigx3_2) * prC) - (2 * sigx14_1 + 2 *
      sigx14_2)^2/sigx7)/sigx7, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((dmusig1 * ewv1_h/(4 *
      (ewu * ewu_h)) - 2 * (ewu * pmusig1/(2 * ewu)^2)) *
      ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) * dmusig1 *
      ewv1_h/ewu_h + (0.5 * sigx9 + 2 * sigx10_1) * sigx16_1)) *
      sigx2_1 - (2 * ((depsi1 * ewv1_h/ewu_h - pepsi1) *
      ewv1/ewu) - ((epsivu1 * sigx17_1 - 0.5) * depsi1 *
      ewv1_h/ewu_h + sigx18_1 * sigx12_1)) * sigx3_1)/(wzdeno *
      sigx7) - wzdeno * (2 * sigx14_1 + 2 * sigx14_2) *
      (sigx2_1 * sigx16_1 - sigx18_1 * sigx3_1)/(wzdeno *
      sigx7)^2) * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * ((((dmusig2 * ewv2_h/(4 *
      (ewu * ewu_h)) - 2 * (ewu * pmusig2/(2 * ewu)^2)) *
      ewv2 - ((0.5 * (sigx15_2 * musig2) - 0.25) * dmusig2 *
      ewv2_h/ewu_h + (0.5 * sigx9 + 2 * sigx10_2) * sigx16_2)) *
      sigx2_2 - ((2 * sigx14_1 + 2 * sigx14_2) * (sigx2_2 *
      sigx16_2 - sigx18_2 * sigx3_2)/sigx7 + (2 * ((depsi2 *
      ewv2_h/ewu_h - pepsi2) * ewv2/ewu) - ((epsivu2 *
      sigx17_2 - 0.5) * depsi2 * ewv2_h/ewu_h + sigx18_2 *
      sigx12_2)) * sigx3_2)) * prC/sigx7), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx11_1 * sigx2_1 -
      sigx13_1 * sigx3_1) - 2 * (sigx11_2 * sigx2_2 - sigx13_2 *
      sigx3_2))/(wzdeno * sigx7) - wzdeno * (2 * sigx14_1 +
      2 * sigx14_2) * (2 * sigx6_1 - 2 * sigx19)/(wzdeno *
      sigx7)^2) * prC * ewz, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((sigx16_1/2 + (pmusig1 -
      sigx15_1 * dmusig1)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx15_1^2 * musig1) *
      dmusig1) * sigx2_1 - ((2 * sigx18_1 + 2 * (pepsi1 -
      depsi1 * sigx17_1)) * ewv1/ewu - (0.25 * (S * (epsilon)/ewv1_h) +
      0.5 * (ewv1_h/ewu_h) - epsivu1 * sigx17_1^2) * depsi1) *
      sigx3_1)/(wzdeno * sigx7) - 2 * ((sigx2_1 * sigx16_1 -
      sigx18_1 * sigx3_1)^2 * ewz/(wzdeno * sigx7)^2)) *
      ewz), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (4 * (prC * (sigx2_1 * sigx16_1 - sigx18_1 * sigx3_1) *
      (sigx2_2 * sigx16_2 - sigx18_2 * sigx3_2) * ewz/(wzdeno *
      sigx7^2))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (2/(wzdeno * sigx7) -
      2 * ((2 * sigx6_1 - 2 * sigx19) * ewz/(wzdeno * sigx7)^2)) *
      (sigx2_1 * sigx16_1 - sigx18_1 * sigx3_1) * ewz,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * 2 * ((((sigx16_2/2 + (pmusig2 - sigx15_2 *
      dmusig2)/2) * ewv2/ewu - (0.25 * (ewv2_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv2_h) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2 - (((2 * sigx18_2 + 2 * (pepsi2 -
      depsi2 * sigx17_2)) * ewv2/ewu - (0.25 * (S * (epsilon)/ewv2_h) +
      0.5 * (ewv2_h/ewu_h) - epsivu2 * sigx17_2^2) * depsi2) *
      sigx3_2 + 2 * (prC * (sigx2_2 * sigx16_2 - sigx18_2 *
      sigx3_2)^2/sigx7))) * prC/sigx7), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (prC * (2 * (prC * wzdeno *
      (2 * sigx6_1 - 2 * sigx19)/(wzdeno * sigx7)^2) +
      2/(wzdeno * sigx7)) * (sigx2_2 * sigx16_2 - sigx18_2 *
      sigx3_2) * ewz), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (prC/(wzdeno * sigx7) - (2 *
      prC + 2 * (ewz/wzdeno)) * sigx6_1 * ewz/(wzdeno *
      sigx7)^2) * prC * (2 * sigx6_1 - 2 * sigx19) * ewz,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# different sigma_u
chessmcesfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  epsi1 <- (2 * (ewv1_h/ewu1_h) + S * (epsilon)/ewv1_h)
  epsi2 <- (2 * (ewv2_h/ewu2_h) + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsi1, 0, 1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  pepsi2 <- pnorm(-epsi2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu1_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu2_h)
  sigx2_1 <- (depsi1/ewv1_h - 2 * (pepsi1/ewu1_h))
  sigx2_2 <- (depsi2/ewv2_h - 2 * (pepsi2/ewu2_h))
  sigx3_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx3_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx4_1 <- exp(2 * (ewv1/ewu1) + 2 * (S * (epsilon)/ewu1_h))
  sigx4_2 <- exp(2 * (ewv2/ewu2) + 2 * (S * (epsilon)/ewu2_h))
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
  sigx10_1 <- (dmusig1 * ewv1_h/ewu1_h)
  sigx10_2 <- (dmusig2 * ewv2_h/ewu2_h)
  sigx11_1 <- (0.5 * (S * (epsilon)/ewu1_h) + 2 * (ewu1 * ewv1/(2 *
    ewu1)^2))
  sigx11_2 <- (0.5 * (S * (epsilon)/ewu2_h) + 2 * (ewu2 * ewv2/(2 *
    ewu2)^2))
  uvsi1 <- (2 * (ewv1/ewu1) + S * (epsilon)/ewu1_h)
  uvsi2 <- (2 * (ewv2/ewu2) + S * (epsilon)/ewu2_h)
  sigx12_1 <- (depsi1 * ewv1_h/ewu1_h - uvsi1 * pepsi1)
  sigx12_2 <- (depsi2 * ewv2_h/ewu2_h - uvsi2 * pepsi2)
  sp1 <- (0.5 * sigx10_1 - sigx11_1 * pmusig1)
  sp2 <- (0.5 * sigx10_2 - sigx11_2 * pmusig2)
  sigx13_1 <- (sp1 * sigx3_1 - sigx12_1 * sigx4_1)
  sigx13_2 <- (sp2 * sigx3_2 - (sigx12_2 * sigx4_2 + 0.5 *
    sigx7_2))
  sigx14_1 <- (wzdeno * ewu1_h * sigx7_1/(wzdeno * ewu1_h)^2)
  sigx15_1 <- (sigx13_1/(wzdeno * ewu1_h) - 0.5 * sigx14_1)
  sigx16_1 <- (0.5 * (ewv1_h/ewu1_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx16_2 <- (0.5 * (ewv2_h/ewu2_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx17_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx16_1 * dmusig1)
  sigx17_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx16_2 * dmusig2)
  sigx18_1 <- (ewv1_h/ewu1_h - 0.5 * (S * (epsilon)/ewv1_h))
  sigx18_2 <- (ewv2_h/ewu2_h - 0.5 * (S * (epsilon)/ewv2_h))
  sigx19_1 <- (sigx3_1 * sigx17_1 - (2 * (ewv1 * pepsi1/ewu1) -
    depsi1 * sigx18_1) * sigx4_1)
  sigx19_2 <- (sigx3_2 * sigx17_2 - (2 * (ewv2 * pepsi2/ewu2) -
    depsi2 * sigx18_2) * sigx4_2)
  sigx20_1 <- (1/(wzdeno * ewu1_h) - ewu1_h * ewz/(wzdeno *
    ewu1_h)^2)
  sigx20_2 <- (prC * sigx7_2/(wzdeno * ewu2_h))
  sigx21 <- (2 * (sigx20_1 * sigx7_1) - 2 * sigx20_2)
  wzsig <- (wzdeno * s8 * ewu1_h)
  sigx22_1 <- (0.5 - wzdeno^2 * ewu1_h^2/(wzdeno * ewu1_h)^2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (2 * ((((musig1/ewv1_h - 1/ewu1_h) *
      dmusig1/ewv1_h - sigx1_1/ewu1_h) * sigx3_1 - ((epsi1/ewv1_h -
      2/ewu1_h) * depsi1/ewv1_h - 2 * (sigx2_1/ewu1_h)) *
      sigx4_1) * ewz/(wzdeno * ewu1_h)) + 2 * ((((musig2/ewv2_h -
      1/ewu2_h) * dmusig2/ewv2_h - sigx1_2/ewu2_h) * sigx3_2 -
      ((epsi2/ewv2_h - 2/ewu2_h) * depsi2/ewv2_h - 2 *
        (sigx2_2/ewu2_h)) * sigx4_2) * prC/ewu2_h) -
      (2 * sigx6_1 + 2 * sigx6_2)^2/s8)/s8, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 *
      (S * (epsilon)/ewu1_h) + 2 * (ewu1 * ewv1/(2 * ewu1)^2)) *
      pmusig1 - 0.5 * sigx10_1)/ewu1_h + (0.5 * (musig1/ewu1_h) -
      sigx11_1/ewv1_h) * dmusig1) * sigx3_1 - ((epsi1/ewu1_h -
      uvsi1/ewv1_h) * depsi1 + (pepsi1 - 2 * sigx12_1)/ewu1_h) *
      sigx4_1)/(wzdeno * ewu1_h) - (sigx15_1 * sigx9 +
      0.5 * (sigx5_1 * wzdeno * ewu1_h/(wzdeno * ewu1_h)^2))) *
      ewz/s8), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 *
      (S * (epsilon)/ewu2_h) + 2 * (ewu2 * ewv2/(2 * ewu2)^2)) *
      pmusig2 - 0.5 * sigx10_2)/ewu2_h + (0.5 * (musig2/ewu2_h) -
      sigx11_2/ewv2_h) * dmusig2) * sigx3_2 - (((epsi2/ewu2_h -
      uvsi2/ewv2_h) * depsi2 + (pepsi2 - 2 * sigx12_2)/ewu2_h) *
      sigx4_2 + 0.5 * sigx5_2))/(s8 * ewu2_h) - sigx13_2 *
      (2 * sigx6_1 + 2 * sigx6_2) * ewu2_h/(s8 * ewu2_h)^2) *
      prC), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((dmusig1 * (ewv1/(2 * ewu1) - (sigx16_1 *
    musig1 + 0.5))/ewv1_h - sigx17_1/ewu1_h) * sigx3_1 -
    ((2 * (ewv1/ewu1) - (epsi1 * sigx18_1 + 0.5)) * depsi1/ewv1_h -
      2 * ((2 * (ewv1 * pepsi1/ewu1) - depsi1 * sigx18_1)/ewu1_h)) *
      sigx4_1)/wzsig - wzdeno * (2 * sigx6_1 + 2 * sigx6_2) *
    ewu1_h * sigx19_1/wzsig^2) * ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((dmusig2 * (ewv2/(2 *
      ewu2) - (sigx16_2 * musig2 + 0.5))/ewv2_h - sigx17_2/ewu2_h) *
      sigx3_2 - ((2 * (ewv2/ewu2) - (epsi2 * sigx18_2 +
      0.5)) * depsi2/ewv2_h - 2 * ((2 * (ewv2 * pepsi2/ewu2) -
      depsi2 * sigx18_2)/ewu2_h)) * sigx4_2)/(s8 * ewu2_h) -
      (2 * sigx6_1 + 2 * sigx6_2) * ewu2_h * sigx19_2/(s8 *
        ewu2_h)^2) * prC), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (sigx5_1 * sigx20_1) -
      ((2 * sigx6_1 + 2 * sigx6_2) * sigx21/s8 + 2 * (sigx5_2 *
        prC/(wzdeno * ewu2_h)))) * ewz/s8, FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((0.5 * (0.5 * (ewv1_h * musig1/ewu1_h) - 0.5) -
    0.5 * sigx11_1) * dmusig1 * ewv1_h/ewu1_h - (sp1 * sigx11_1 +
    (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) * ewu1 * ewv1/(2 *
      ewu1)^2) - 0.25 * (S * (epsilon)/ewu1_h)) * pmusig1)) *
    sigx3_1 - (((epsi1 * ewv1_h - S * (epsilon))/ewu1_h -
    (0.5 + 2 * (ewv1/ewu1))) * depsi1 * ewv1_h/ewu1_h + (0.5 *
    (S * (epsilon)/ewu1_h) + 2 * (ewv1/ewu1)) * pepsi1 -
    uvsi1 * sigx12_1) * sigx4_1)/(wzdeno * ewu1_h) - ((0.5 *
    (sigx22_1 * sigx7_1 + sp1 * sigx3_1 - sigx12_1 * sigx4_1) +
    0.5 * sigx13_1) * wzdeno * ewu1_h/(wzdeno * ewu1_h)^2 +
    2 * (sigx15_1^2 * ewz/s8))) * ewz/s8), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (sigx15_1 * sigx13_2 * prC * ewu2_h *
      ewz/(s8 * ewu2_h)^2)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((dmusig1 * ewv1_h/(4 *
      (ewu1 * ewu1_h)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (sigx16_1 * musig1) - 0.25) * dmusig1 *
      ewv1_h/ewu1_h + sigx11_1 * sigx17_1)) * sigx3_1 -
      (2 * ((depsi1 * ewv1_h/ewu1_h - pepsi1) * ewv1/ewu1) -
        ((epsi1 * sigx18_1 - 0.5) * depsi1 * ewv1_h/ewu1_h +
          (2 * (ewv1 * pepsi1/ewu1) - depsi1 * sigx18_1) *
          uvsi1)) * sigx4_1)/wzsig - (0.5 * s8 + 2 *
      (sigx15_1 * ewz)) * wzdeno * ewu1_h * sigx19_1/wzsig^2) *
      ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (sigx15_1 * prC * ewu2_h *
      sigx19_2 * ewz/(s8 * ewu2_h)^2)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * (sigx13_1 * sigx20_1 - (sigx22_1 * ewz + 0.5 * wzdeno) *
      ewu1_h * sigx7_1/(wzdeno * ewu1_h)^2) - 2 * (sigx15_1 *
      sigx21 * ewz/s8)) * ewz/s8, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((0.5 * (0.5 * (ewv2_h *
      musig2/ewu2_h) - 0.5) - 0.5 * sigx11_2) * dmusig2 *
      ewv2_h/ewu2_h - (sp2 * sigx11_2 + (2 * ((1 - 8 *
      (ewu2^2/(2 * ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) -
      0.25 * (S * (epsilon)/ewu2_h)) * pmusig2)) * sigx3_2 -
      ((((epsi2 * ewv2_h - S * (epsilon))/ewu2_h - (0.5 +
        2 * (ewv2/ewu2))) * depsi2 * ewv2_h/ewu2_h +
        (0.5 * (S * (epsilon)/ewu2_h) + 2 * (ewv2/ewu2)) *
          pepsi2 - uvsi2 * sigx12_2) * sigx4_2 + 0.5 *
        (sp2 * sigx3_2 - sigx12_2 * sigx4_2)))/(s8 *
      ewu2_h) - sigx13_2 * (0.5 * (s8 * ewu2_h) + 2 * (sigx13_2 *
      prC))/(s8 * ewu2_h)^2) * prC), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (sigx13_2 * prC * wzdeno *
      ewu1_h * sigx19_1 * ewz/(wzsig^2 * ewu2_h))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((dmusig2 * ewv2_h/(4 * (ewu2 * ewu2_h)) - 2 *
    (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (sigx16_2 *
    musig2) - 0.25) * dmusig2 * ewv2_h/ewu2_h + sigx11_2 *
    sigx17_2)) * sigx3_2 - (2 * ((depsi2 * ewv2_h/ewu2_h -
    pepsi2) * ewv2/ewu2) - ((epsi2 * sigx18_2 - 0.5) * depsi2 *
    ewv2_h/ewu2_h + (2 * (ewv2 * pepsi2/ewu2) - depsi2 *
    sigx18_2) * uvsi2)) * sigx4_2)/(s8 * ewu2_h) - (0.5 *
    (s8 * ewu2_h) + 2 * (sigx13_2 * prC)) * sigx19_2/(s8 *
    ewu2_h)^2) * prC), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (prC * (2 * (sigx13_2 * sigx21/(s8 *
      ewu2_h)) + 2 * ((sp2 * sigx3_2 - sigx12_2 * sigx4_2)/(wzdeno *
      ewu2_h) - 0.5 * (wzdeno * ewu2_h * sigx7_2/(wzdeno *
      ewu2_h)^2))) * ewz/s8), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((sigx17_1/2 + (pmusig1 -
      sigx16_1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewv1_h/ewu1_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx16_1^2 * musig1) *
      dmusig1) * sigx3_1 - ((2 * (2 * (ewv1 * pepsi1/ewu1) -
      depsi1 * sigx18_1) + 2 * (pepsi1 - depsi1 * sigx18_1)) *
      ewv1/ewu1 - (0.25 * (S * (epsilon)/ewv1_h) + 0.5 *
      (ewv1_h/ewu1_h) - epsi1 * sigx18_1^2) * depsi1) *
      sigx4_1)/wzsig - 2 * (sigx19_1^2 * ewz/wzsig^2)) *
      ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (prC * ewu2_h * sigx19_1 * sigx19_2 *
      ewz/((s8 * ewu2_h)^2 * wzdeno * ewu1_h))), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * sigx20_1 - 2 * (sigx21 *
      ewz/wzsig)) * sigx19_1 * ewz/s8, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((sigx17_2/2 + (pmusig2 -
      sigx16_2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewv2_h/ewu2_h) +
      0.25 * (S * (epsilon)/ewv2_h) - sigx16_2^2 * musig2) *
      dmusig2) * sigx3_2 - ((2 * (2 * (ewv2 * pepsi2/ewu2) -
      depsi2 * sigx18_2) + 2 * (pepsi2 - depsi2 * sigx18_2)) *
      ewv2/ewu2 - (0.25 * (S * (epsilon)/ewv2_h) + 0.5 *
      (ewv2_h/ewu2_h) - epsi2 * sigx18_2^2) * depsi2) *
      sigx4_2)/(s8 * ewu2_h) - 2 * (prC * sigx19_2^2/(s8 *
      ewu2_h)^2)) * prC), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (prC * (2 * (sigx21/s8) +
      2/wzdeno) * sigx19_2 * ewz/(s8 * ewu2_h)), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    ((2 * (prC * (1/(wzdeno^2 * ewu2_h) + ewu2_h/(wzdeno *
      ewu2_h)^2) * sigx7_2) - (sigx21^2/s8 + 2 * ((2 -
      2 * (wzdeno * ewu1_h^2 * ewz/(wzdeno * ewu1_h)^2)) *
      ewu1_h * sigx7_1/(wzdeno * ewu1_h)^2))) * ewz + 2 *
      (sigx20_1 * sigx7_1) - 2 * sigx20_2) * ewz/s8, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf generalized exponential-normal distribution
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
cnsfgenexponormAlgOpt <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfgenexponormlike(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(ccnsfgenexponormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfgenexponormlike, grad = cgradcnsfgenexponormlike,
      hess = chesscnsfgenexponormlike, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfgenexponormlike(mleObj$par,
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
      mleObj$hessian <- chesscnsfgenexponormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfgenexponormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfgenexponormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfgenexponormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGenExpo = initGenExpo))
}

# different sigma_u
mcesfgenexponormAlgOpt <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfgenexponormlike(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmcesfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfgenexponormlike, grad = cgradmcesfgenexponormlike,
      hess = chessmcesfgenexponormlike, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfgenexponormlike(mleObj$par,
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
      mleObj$hessian <- chessmcesfgenexponormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfgenexponormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfgenexponormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfgenexponormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGenExpo = initGenExpo))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf genexpo-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u
ccnsfgenexponormeff <- function(object, level) {
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
  A1 <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B1 <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a1 <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b1 <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  A2 <- object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))
  B2 <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv2)/exp(Wu)
  a2 <- -object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)
  b2 <- -object$S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(exp(A1) * pnorm(a1) -
    exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(exp(A2) * pnorm(a2) -
    exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
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
cmcesfgenexponormeff <- function(object, level) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv1)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv2)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(exp(A1) * pnorm(a1) -
    exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(exp(A2) * pnorm(a2) -
    exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
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
#' marginal impact on efficiencies for cnsf genexpo-normal distribution
#' @param object object of class sfacross
#' @noRd
# same sigma_u
ccnsfmarggenexponorm_Eu <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B1 <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a1 <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b1 <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  A2 <- object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))
  B2 <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv2)/exp(Wu)
  a2 <- -object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)
  b2 <- -object$S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmarggenexponorm_Vu <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B1 <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a1 <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b1 <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  A2 <- object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))
  B2 <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv2)/exp(Wu)
  a2 <- -object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)
  b2 <- -object$S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# different sigma_u
cmcesfmarggenexponorm_Eu <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv1)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv2)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmarggenexponorm_Vu <- function(object) {
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
  A1 <- object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))
  B1 <- 2 * object$S * epsilon/exp(Wu1/2) + 2 * exp(Wv1)/exp(Wu1)
  a1 <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)
  b1 <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu1/2)
  A2 <- object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))
  B2 <- 2 * object$S * epsilon/exp(Wu2/2) + 2 * exp(Wv2)/exp(Wu2)
  a2 <- -object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)
  b2 <- -object$S * epsilon/exp(Wv2/2) - 2 * exp(Wv2/2)/exp(Wu2/2)
  Pi1 <- 2/exp(Wu1/2) * (exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- 2/exp(Wu2/2) * (exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
