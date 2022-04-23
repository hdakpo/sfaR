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
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf halfnormal-normal distribution
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
ccnsfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar1 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# different sigma_u
cmcesfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf halfnormal-normal distribution
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
cstcnsfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + halfnormal - normal distributions...\n")
  initHalf <- maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradhalfnormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initHalf$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("CN_", colnames(Zvar)))
  names(initHalf$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# different sigma_u
cstmcesfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("nitialization: SFA + halfnormal - normal distributions...\n")
  initHalf <- maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradhalfnormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = printInfo, reltol = tol), nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]),
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar)
  Esti <- initHalf$estimate
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("MCE_",
    colnames(Zvar)))
  names(initHalf$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf halfnormal-normal distribution
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
cgradcnsfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  depsi1 <- dnorm(S * (epsilon)/sqrt(sigma_sq1))
  depsi2 <- dnorm(S * (epsilon)/sqrt(sigma_sq2))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewvsq1 <- (1 - ewv1/(sigma_sq1))
  ewvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx1_1 <- (0.5 * (ewvsq1 * ewu/sigmastar1) + sigmastar1)
  sigx1_2 <- (0.5 * (ewvsq2 * ewu/sigmastar2) + sigmastar2)
  sigx2_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx2_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx3_1 <- (sigx1_1 * dmusig1 * depsi1 * ewu/ssq1^2 + 0.5 *
    sigx2_1)
  sigx3_2 <- (sigx1_2 * dmusig2 * depsi2 * ewu/ssq2^2 + 0.5 *
    sigx2_2)
  sigx4_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon))
  sigx4_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + S * depsi2 *
    pmusig2 * (epsilon))
  sisi1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sisi2 <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx5_1 <- (sigx4_1 * ewz/sisi1)
  sigx5_2 <- (prC * sigx4_2/sisi2)
  wzsq1 <- (wzdeno * sqrt(sigma_sq1))
  wzsq2 <- (wzdeno * sqrt(sigma_sq2))
  sigx6_1 <- (depsi1 * ewz * pmusig1/wzsq1)
  sigx6_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  prU1 <- (1 - ewu/(sigma_sq1))
  prU2 <- (1 - ewu/(sigma_sq2))
  sigx7_1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (wzdeno * depsi1 * pmusig1/wzsq1^2)
  sigx8_2 <- (prC * depsi2 * pmusig2/wzsq2)
  sigx9_1 <- (1/wzsq1 - ewz * sqrt(sigma_sq1)/wzsq1^2)
  sigx9_2 <- ((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq2))
  sigx10_1 <- (sigx9_1 * depsi1 * pmusig1)
  sigx10_2 <- (depsi2 * pmusig2/(sigma_sq2))
  ssqx1 <- (1/ssq1 - sigx7_1 * ewu/ssq1^2)
  ssqx2 <- (1/ssq2 - sigx7_2 * ewu/ssq2^2)
  sigx11_1 <- (0.5 * sigx2_1 - ssqx1 * dmusig1 * depsi1)
  sigx11_2 <- (0.5 * sigx2_2 - ssqx2 * dmusig2 * depsi2)
  sigx12_1 <- (S * sigx11_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx12_2 <- (S * sigx11_2 * (epsilon) - 0.5 * sigx10_2)
  sigx13_1 <- (S * sigx3_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx13_2 <- (S * sigx3_2 * (epsilon) - 0.5 * sigx10_2)
  sigx14_1 <- (ewz * sigx12_1/sqrt(sigma_sq1))
  sigx14_2 <- (prC * sigx12_2/sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (2 *
    sigx5_2 + 2 * sigx5_1)/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (2 * sigx14_2 + 2 *
      sigx14_1) * ewu/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (ewv1 * ewz * sigx13_1/((2 *
      sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1))), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (prC * ewv2 * sigx13_2/sigx9_2),
      FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 *
      sigx10_1 - 2 * sigx8_2) * ewz/(2 * sigx6_2 + 2 *
      sigx6_1), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u
cgradmcesfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  epsi1 <- S * (epsilon)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon)/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sxq <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 *
    pmusig2 * (epsilon))
  wzsq1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  wzxsq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (2 * (prC * sigx1_2/sxq) + 2 * (sigx1_1 * ewz/wzsq1))
  sigx3_1 <- (depsi1 * ewz * pmusig1/wzxsq1)
  sigx3_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx4_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx4_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  usq1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  usq2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
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
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  sigx9_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  s9sq2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx9_2 <- (s9sq2 * dmusig2 * depsi2 * ewu2/ssq2^2 + 0.5 *
    sigx4_2)
  sigx10_1 <- (sigx9_1 * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 *
    sigx4_1)
  sigx10_2 <- (S * sigx9_2 * (epsilon) - 0.5 * s7sq2)
  s12sq1 <- (1/wzxsq1 - ewz * sqrt(sigma_sq1)/wzxsq1^2)
  sigx12 <- (2 * (s12sq1 * depsi1 * pmusig1) - 2 * sigx11_2)
  sigx18_1 <- (S * sigx10_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/(2 *
    sigx3_2 + 2 * sigx3_1), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = 2 * (ewu1 * ewz * sigx7_1/sigx8_1), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (prC * ewu2 * sigx7_2/sigx8_2),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
      (ewv1 * ewz * sigx18_1/sigx8_1), FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = 2 * (prC * ewv2 * sigx10_2/sigx8_2),
      FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx12 *
      ewz/(2 * sigx3_2 + 2 * sigx3_1), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf halfnormal-normal distribution
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
chesscnsfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  depsi1 <- dnorm(S * (epsilon)/sqrt(sigma_sq1))
  depsi2 <- dnorm(S * (epsilon)/sqrt(sigma_sq2))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewvsq1 <- (1 - ewv1/(sigma_sq1))
  ewvsq2 <- (1 - ewv2/(sigma_sq2))
  sigx1_1 <- (0.5 * (ewvsq1 * ewu/sigmastar1) + sigmastar1)
  sigx1_2 <- (0.5 * (ewvsq2 * ewu/sigmastar2) + sigmastar2)
  sigx2_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx2_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx3_1 <- (sigx1_1 * dmusig1 * depsi1 * ewu/ssq1^2 + 0.5 *
    sigx2_1)
  sigx3_2 <- (sigx1_2 * dmusig2 * depsi2 * ewu/ssq2^2 + 0.5 *
    sigx2_2)
  sigx4_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon))
  sigx4_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + S * depsi2 *
    pmusig2 * (epsilon))
  sisi1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sisi2 <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx5_1 <- (sigx4_1 * ewz/sisi1)
  sigx5_2 <- (prC * sigx4_2/sisi2)
  wzsq1 <- (wzdeno * sqrt(sigma_sq1))
  wzsq2 <- (wzdeno * sqrt(sigma_sq2))
  sigx6_1 <- (depsi1 * ewz * pmusig1/wzsq1)
  sigx6_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  prU1 <- (1 - ewu/(sigma_sq1))
  prU2 <- (1 - ewu/(sigma_sq2))
  sigx7_1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (wzdeno * depsi1 * pmusig1/wzsq1^2)
  sigx8_2 <- (prC * depsi2 * pmusig2/wzsq2)
  sigx9_1 <- (1/wzsq1 - ewz * sqrt(sigma_sq1)/wzsq1^2)
  sigx9_2 <- ((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq2))
  sigx10_1 <- (sigx9_1 * depsi1 * pmusig1)
  sigx10_2 <- (depsi2 * pmusig2/(sigma_sq2))
  ssqx1 <- (1/ssq1 - sigx7_1 * ewu/ssq1^2)
  ssqx2 <- (1/ssq2 - sigx7_2 * ewu/ssq2^2)
  sigx11_1 <- (0.5 * sigx2_1 - ssqx1 * dmusig1 * depsi1)
  sigx11_2 <- (0.5 * sigx2_2 - ssqx2 * dmusig2 * depsi2)
  sigx12_1 <- (S * sigx11_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx12_2 <- (S * sigx11_2 * (epsilon) - 0.5 * sigx10_2)
  sigx13_1 <- (S * sigx3_1 * (epsilon)/wzdeno - 0.5 * sigx8_1)
  sigx13_2 <- (S * sigx3_2 * (epsilon) - 0.5 * sigx10_2)
  sigx14_1 <- (ewz * sigx12_1/sqrt(sigma_sq1))
  sigx14_2 <- (prC * sigx12_2/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * ewu/sigmastar1 + S * pmusig1 * (epsilon))
  sigx15_2 <- (dmusig2 * ewu/sigmastar2 + S * pmusig2 * (epsilon))
  ssqxx1 <- (S * sigx15_1 * (epsilon)/(sigma_sq1) - pmusig1)
  ssqxx2 <- (S * sigx15_2 * (epsilon)/(sigma_sq2) - pmusig2)
  dd1 <- (depsi1 * ewu/ewv1 + depsi1)
  dd2 <- (depsi2 * ewu/ewv2 + depsi2)
  sigx17_1 <- (0.5 * ssqxx1 - 0.5 * pmusig1)
  sigx17_2 <- (0.5 * ssqxx2 - 0.5 * pmusig2)
  sigx18_1 <- (wzdeno * sigx4_1/(wzsq1^2 * (sigma_sq1)))
  sigx18_2 <- (sigx4_2/(sigma_sq2))
  pesig1 <- (S * pmusig1 * (epsilon)/(sigma_sq1)^2)
  pesig2 <- (S * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx19_1 <- (0.5 * pesig1 - ssqx1 * dmusig1)
  sigx19_2 <- (0.5 * pesig2 - ssqx2 * dmusig2)
  sigx20_1 <- (S * sigx19_1 * (epsilon) - 2 * (pmusig1/(sigma_sq1)))
  sigx20_2 <- (S * sigx19_2 * (epsilon) - 2 * (pmusig2/(sigma_sq2)))
  sigx21_1 <- (S * depsi1 * sigx20_1 * (epsilon)/(sigma_sq1)^2)
  sigx21_2 <- (S * depsi2 * sigx20_2 * (epsilon)/(sigma_sq2)^2)
  sigx22_1 <- (S * sigx11_1 * (epsilon) - wzdeno^2 * depsi1 *
    pmusig1/wzsq1^2)
  sigx22_2 <- (S * sigx11_2 * (epsilon) - depsi2 * pmusig2/(sigma_sq2))
  sigx23_1 <- (wzdeno * sigx22_1/wzsq1^2)
  sigx23_2 <- (wzdeno * depsi2 * pmusig2/wzsq2^2)
  sigx24_1 <- (0.5/sqrt(sigma_sq1) - wzdeno^2 * sqrt(sigma_sq1)/wzsq1^2)
  sigx25 <- (2 * sigx10_1 - 2 * sigx8_2)/(2 * sigx6_2 + 2 *
    sigx6_1)
  sigx26_1 <- (0.5 * (S^2 * sigx9_1 * depsi1 * (epsilon)^2/(sigma_sq1)^2) -
    (sigx24_1 * ewz + 0.5 * (wzdeno/sqrt(sigma_sq1))) * depsi1/wzsq1^2)
  sigx27 <- ((2 * sigx6_2 + 2 * sigx6_1)/sqrt(sigma_sq1))
  sigx28 <- ((2 * sigx6_2 + 2 * sigx6_1)/sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (2 * (prC * (depsi2 * ssqxx2 +
      S * dmusig2 * dd2 * ewu * (epsilon)/ssq2)/((sigma_sq2) *
      sqrt(sigma_sq2))) + 2 * ((depsi1 * ssqxx1 + S * dmusig1 *
      dd1 * ewu * (epsilon)/ssq1) * ewz/(wzdeno * (sigma_sq1) *
      sqrt(sigma_sq1))) - (2 * sigx5_2 + 2 * sigx5_1)^2/(2 *
      sigx6_2 + 2 * sigx6_1))/(2 * sigx6_2 + 2 * sigx6_1),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (((ssqx1 * dmusig1 *
      depsi1 + S * (sigx17_1 * depsi1/(sigma_sq1) - S *
      ssqx1 * dmusig1 * dd1 * (epsilon)) * (epsilon)/(sigma_sq1))/wzdeno -
      0.5 * sigx18_1) * ewz/sqrt(sigma_sq1)) + 2 * ((ssqx2 *
      dmusig2 * depsi2 + (S * (sigx17_2 * depsi2/(sigma_sq2) -
      S * ssqx2 * dmusig2 * dd2 * (epsilon)) * (epsilon) -
      0.5 * sigx18_2)/(sigma_sq2)) * prC/sqrt(sigma_sq2)) -
      (2 * sigx5_2 + 2 * sigx5_1) * (2 * sigx14_2 + 2 *
        sigx14_1)/(2 * sigx6_2 + 2 * sigx6_1)) * ewu/(2 *
      sigx6_2 + 2 * sigx6_1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((S * (sigx17_1 * depsi1/(sigma_sq1) + S *
    sigx1_1 * dmusig1 * dd1 * ewu * (epsilon)/ssq1^2) * (epsilon)/(sigma_sq1) -
    sigx1_1 * dmusig1 * depsi1 * ewu/ssq1^2)/wzdeno - 0.5 *
    sigx18_1)/((2 * sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1)) -
    (2 * sigx5_2 + 2 * sigx5_1) * sigx13_1 * sqrt(sigma_sq1)/((2 *
      sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1))^2) * ewv1 *
    ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * (((S * (sigx17_2 * depsi2/(sigma_sq2) +
      S * sigx1_2 * dmusig2 * dd2 * ewu * (epsilon)/ssq2^2) *
      (epsilon) - 0.5 * sigx18_2)/(sigma_sq2) - sigx1_2 *
      dmusig2 * depsi2 * ewu/ssq2^2)/sigx9_2 - (2 * sigx5_2 +
      2 * sigx5_1) * sigx13_2 * sqrt(sigma_sq2)/sigx9_2^2) *
      prC * ewv2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (sigx9_1 * sigx4_1/(sigma_sq1)) -
      ((2 * sigx5_2 + 2 * sigx5_1) * sigx25 + 2 * (prC *
        sigx4_2/(wzdeno * (sigma_sq2) * sqrt(sigma_sq2))))) *
      ewz/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (prC * (S * (0.5 * sigx21_2 - (0.5 * (S^2 * ssqx2 *
      depsi2 * (epsilon)^2/(sigma_sq2)^2) - (((0.5 * (ewu/(sigma_sq2)) +
      1 - 0.5 * (0.5 * prU2 + ewu/(sigma_sq2))) * prU2 *
      ewv2/sigmastar2 + (2 - 2 * (sigx7_2^2 * ewu * (sigma_sq2)/ssq2^2)) *
      sigmastar2)/ssq2^2 + S^2 * ssqx2^2 * ewu * (epsilon)^2/ssq2) *
      depsi2) * dmusig2) * (epsilon) - (0.5 * sigx12_2 +
      0.5 * sigx22_2)/(sigma_sq2))/sqrt(sigma_sq2)) + 2 *
      (ewz * (S * (0.5 * sigx21_1 - (0.5 * (S^2 * ssqx1 *
        depsi1 * (epsilon)^2/(sigma_sq1)^2) - (((0.5 *
        (ewu/(sigma_sq1)) + 1 - 0.5 * (0.5 * prU1 + ewu/(sigma_sq1))) *
        prU1 * ewv1/sigmastar1 + (2 - 2 * (sigx7_1^2 *
        ewu * (sigma_sq1)/ssq1^2)) * sigmastar1)/ssq1^2 +
        S^2 * ssqx1^2 * ewu * (epsilon)^2/ssq1) * depsi1) *
        dmusig1) * (epsilon)/wzdeno - (0.5 * sigx23_1 +
        0.5 * (sigx12_1/(sigma_sq1))))/sqrt(sigma_sq1)) -
      (2 * sigx14_2 + 2 * sigx14_1)^2/(2 * sigx6_2 + 2 *
        sigx6_1)) * ewu + 2 * sigx14_2 + 2 * sigx14_1) *
    ewu/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU1 *
      ewv1) - S^2 * sigx1_1 * ssqx1 * ewu * (epsilon)^2)/(sigma_sq1) +
      0.5 * ((ewu/(sigma_sq1) - 1) * ewv1/(sigma_sq1) +
        1 - 0.5 * (prU1 * ewvsq1))) * depsi1/sigmastar1 +
      0.5 * (S^2 * sigx1_1 * depsi1 * (epsilon)^2/(sigma_sq1)^2)) *
      ewu + sigx1_1 * (1 - 2 * (sigx7_1 * ewu * (sigma_sq1) *
      sigmastar1/ssq1^2)) * depsi1) * dmusig1/ssq1^2 +
      0.5 * sigx21_1) * (epsilon)/wzdeno - 0.5 * sigx23_1)/((2 *
      sigx6_2 + 2 * sigx6_1) * sqrt(sigma_sq1)) - ((2 *
      sigx14_2 + 2 * sigx14_1) * sqrt(sigma_sq1) + 0.5 *
      sigx27) * sigx13_1/((2 * sigx6_2 + 2 * sigx6_1) *
      sqrt(sigma_sq1))^2) * ewu * ewv1 * ewz), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU2 *
      ewv2) - S^2 * sigx1_2 * ssqx2 * ewu * (epsilon)^2)/(sigma_sq2) +
      0.5 * ((ewu/(sigma_sq2) - 1) * ewv2/(sigma_sq2) +
        1 - 0.5 * (prU2 * ewvsq2))) * depsi2/sigmastar2 +
      0.5 * (S^2 * sigx1_2 * depsi2 * (epsilon)^2/(sigma_sq2)^2)) *
      ewu + sigx1_2 * (1 - 2 * (sigx7_2 * ewu * (sigma_sq2) *
      sigmastar2/ssq2^2)) * depsi2) * dmusig2/ssq2^2 +
      0.5 * sigx21_2) * (epsilon) - 0.5 * (sigx22_2/(sigma_sq2)))/sigx9_2 -
      ((2 * sigx14_2 + 2 * sigx14_1) * sqrt(sigma_sq2) +
        0.5 * sigx28) * sigx13_2/sigx9_2^2) * prC * ewu *
      ewv2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx26_1 * pmusig1 -
      S * sigx9_1 * ssqx1 * dmusig1 * depsi1 * (epsilon)) -
      ((2 * sigx14_2 + 2 * sigx14_1) * sigx25 + 2 * (prC *
        (S * sigx11_2 * (epsilon)/wzdeno - 0.5 * sigx23_2)/sqrt(sigma_sq2)))) *
      ewu * ewz/(2 * sigx6_2 + 2 * sigx6_1), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) -
      0.5 * (0.5 * ewvsq1 + ewv1/(sigma_sq1))) * ewvsq1 +
      S^2 * sigx1_1^2 * ewu * ewv1 * (epsilon)^2/(ssq1^2 *
        (sigma_sq1))) * depsi1 * ewu/sigmastar1 + ((0.5 *
      (S^2 * depsi1 * (epsilon)^2/(sigma_sq1)^2) - 2 *
      (sigx1_1 * depsi1 * (sigma_sq1) * sigmastar1/ssq1^2)) *
      ewv1 + depsi1) * sigx1_1) * dmusig1 * ewu/ssq1^2 +
      S * (0.5 * (ewv1 * (S * (sigx1_1 * dmusig1 * ewu/ssq1^2 +
        0.5 * pesig1) * (epsilon) - 2 * (pmusig1/(sigma_sq1)))) +
        0.5 * pmusig1) * depsi1 * (epsilon)/(sigma_sq1)^2) *
      (epsilon)/wzdeno - (0.5 * (depsi1 * pmusig1) + 0.5 *
      (ewv1 * (S * sigx3_1 * (epsilon) - wzdeno^2 * depsi1 *
        pmusig1/wzsq1^2))) * wzdeno/wzsq1^2)/((2 * sigx6_2 +
      2 * sigx6_1) * sqrt(sigma_sq1)) - (0.5 * sigx27 +
      2 * (ewz * sigx13_1)) * ewv1 * sigx13_1/((2 * sigx6_2 +
      2 * sigx6_1) * sqrt(sigma_sq1))^2) * ewv1 * ewz),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (4 * (prC * ewv1 * ewv2 * ewz * sigx13_1 * sigx13_2 *
      sqrt(sigma_sq2)/(sigx9_2^2 * sqrt(sigma_sq1)))),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx26_1 * pmusig1 +
      S * sigx1_1 * sigx9_1 * dmusig1 * depsi1 * ewu *
        (epsilon)/ssq1^2) - 2 * ((2 * sigx10_1 - 2 *
      sigx8_2) * ewz * sigx13_1/((2 * sigx6_2 + 2 * sigx6_1) *
      sqrt(sigma_sq1)))) * ewv1 * ewz/(2 * sigx6_2 + 2 *
      sigx6_1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * 2 * (((S * ((((0.5 * (ewv2/(sigma_sq2)) -
      0.5 * (0.5 * ewvsq2 + ewv2/(sigma_sq2))) * ewvsq2 +
      S^2 * sigx1_2^2 * ewu * ewv2 * (epsilon)^2/(ssq2^2 *
        (sigma_sq2))) * depsi2 * ewu/sigmastar2 + ((0.5 *
      (S^2 * depsi2 * (epsilon)^2/(sigma_sq2)^2) - 2 *
      (sigx1_2 * depsi2 * (sigma_sq2) * sigmastar2/ssq2^2)) *
      ewv2 + depsi2) * sigx1_2) * dmusig2 * ewu/ssq2^2 +
      S * (0.5 * (ewv2 * (S * (sigx1_2 * dmusig2 * ewu/ssq2^2 +
        0.5 * pesig2) * (epsilon) - 2 * (pmusig2/(sigma_sq2)))) +
        0.5 * pmusig2) * depsi2 * (epsilon)/(sigma_sq2)^2) *
      (epsilon) - (0.5 * (depsi2 * pmusig2) + 0.5 * (ewv2 *
      (S * sigx3_2 * (epsilon) - depsi2 * pmusig2/(sigma_sq2))))/(sigma_sq2))/sigx9_2 -
      (0.5 * sigx28 + 2 * (prC * sigx13_2)) * ewv2 * sigx13_2/sigx9_2^2) *
      prC * ewv2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (prC * (2 * ((2 * sigx10_1 -
      2 * sigx8_2) * sigx13_2/(2 * sigx6_2 + 2 * sigx6_1)) +
      2 * (S * sigx3_2 * (epsilon)/wzdeno - 0.5 * sigx23_2)) *
      ewv2 * ewz/sigx9_2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((2 * (prC * (1/(wzdeno^2 *
      sqrt(sigma_sq2)) + sqrt(sigma_sq2)/wzsq2^2) * depsi2 *
      pmusig2) - ((2 * sigx10_1 - 2 * sigx8_2)^2/(2 * sigx6_2 +
      2 * sigx6_1) + 2 * ((2 - 2 * (wzdeno * (sigma_sq1) *
      ewz/wzsq1^2)) * depsi1 * pmusig1 * sqrt(sigma_sq1)/wzsq1^2))) *
      ewz + 2 * sigx10_1 - 2 * sigx8_2) * ewz/(2 * sigx6_2 +
      2 * sigx6_1), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# different sigma_u
chessmcesfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  epsi1 <- S * (epsilon)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon)/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sxq <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 *
    pmusig2 * (epsilon))
  wzsq1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  wzxsq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (2 * (prC * sigx1_2/sxq) + 2 * (sigx1_1 * ewz/wzsq1))
  sigx3_1 <- (depsi1 * ewz * pmusig1/wzxsq1)
  sigx3_2 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx4_1 <- (S * depsi1 * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx4_2 <- (S * depsi2 * pmusig2 * (epsilon)/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  usq1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  usq2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
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
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  sigx9_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  s9sq2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx9_2 <- (s9sq2 * dmusig2 * depsi2 * ewu2/ssq2^2 + 0.5 *
    sigx4_2)
  sigx10_1 <- (sigx9_1 * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 *
    sigx4_1)
  sigx10_2 <- (S * sigx9_2 * (epsilon) - 0.5 * s7sq2)
  s12sq1 <- (1/wzxsq1 - ewz * sqrt(sigma_sq1)/wzxsq1^2)
  sigx12 <- (2 * (s12sq1 * depsi1 * pmusig1) - 2 * sigx11_2)
  sigx18_1 <- (S * sigx10_1 * (epsilon)/wzdeno - 0.5 * sigx11_1)
  sigx13_1 <- (dmusig1 * ewu1/sigmastar1 + S * pmusig1 * (epsilon))
  sigx13_2 <- (dmusig2 * ewu2/sigmastar2 + S * pmusig2 * (epsilon))
  sigx14_1 <- (S * sigx13_1 * (epsilon)/(sigma_sq1) - pmusig1)
  sigx14_2 <- (S * sigx13_2 * (epsilon)/(sigma_sq2) - pmusig2)
  sigx15_1 <- (0.5 * sigx14_1 - 0.5 * pmusig1) * depsi1/(sigma_sq1)
  sigx15_2 <- (0.5 * sigx14_2 - 0.5 * pmusig2) * depsi2/(sigma_sq2)
  sigx16_1 <- (depsi1 * ewu1/ewv1 + depsi1)
  sigx16_2 <- (depsi2 * ewu2/ewv2 + depsi2)
  sigx17_1 <- (wzdeno * sigx1_1/(wzxsq1^2 * (sigma_sq1)))
  sigx17_2 <- (sigx1_2/(sigma_sq2))
  sigx19_1 <- (S^2 * s12sq1 * depsi1 * (epsilon)^2/(sigma_sq1)^2)
  sigx20_1 <- ((0.5/sqrt(sigma_sq1) - wzdeno^2 * sqrt(sigma_sq1)/wzxsq1^2) *
    ewz + 0.5 * (wzdeno/sqrt(sigma_sq1)))
  sigx21_1 <- (0.5 * sigx19_1 - sigx20_1 * depsi1/wzxsq1^2) *
    pmusig1
  sigx22_1 <- (S * pmusig1 * (epsilon)/(sigma_sq1)^2)
  sigx22_2 <- (S * pmusig2 * (epsilon)/(sigma_sq2)^2)
  sigx23_1 <- (S * (0.5 * sigx22_1 - sigx5_1 * dmusig1) * (epsilon) -
    2 * (pmusig1/(sigma_sq1)))
  sigx23_2 <- (S * (0.5 * sigx22_2 - sigx5_2 * dmusig2) * (epsilon) -
    2 * (pmusig2/(sigma_sq2)))
  sigx24_1 <- (S * depsi1 * sigx23_1 * (epsilon)/(sigma_sq1)^2)
  sigx24_2 <- (S * depsi2 * sigx23_2 * (epsilon)/(sigma_sq2)^2)
  sigx25_1 <- (wzdeno * (S * sigx6_1 * (epsilon) - wzdeno^2 *
    depsi1 * pmusig1/wzxsq1^2)/wzxsq1^2)
  sigx25_2 <- ((S * sigx6_2 * (epsilon) - depsi2 * pmusig2/(sigma_sq2))/(sigma_sq2))
  sigx26_1 <- ((2 * sigx3_2 + 2 * sigx3_1)/sqrt(sigma_sq1))
  sigx26_2 <- ((2 * sigx3_2 + 2 * sigx3_1)/sqrt(sigma_sq2))
  sigx27_1 <- (sigx8_2^2 * sqrt(sigma_sq1))
  sigx27_2 <- (wzdeno * depsi2 * pmusig2/(wzdeno * sqrt(sigma_sq2))^2)
  sigx28_1 <- (0.5 * sigx26_1 + 2 * (ewz * sigx7_1))
  sigx28_2 <- (0.5 * sigx26_2 + 2 * (prC * sigx7_2))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (2 * (prC * (depsi2 * sigx14_2 +
      S * dmusig2 * sigx16_2 * ewu2 * (epsilon)/ssq2)/sxq) +
      2 * ((depsi1 * sigx14_1 + S * dmusig1 * sigx16_1 *
        ewu1 * (epsilon)/ssq1) * ewz/wzsq1) - sigx2^2/(2 *
      sigx3_2 + 2 * sigx3_1))/(2 * sigx3_2 + 2 * sigx3_1),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((sigx5_1 * dmusig1 *
      depsi1 + S * (sigx15_1 - S * sigx5_1 * dmusig1 *
      sigx16_1 * (epsilon)) * (epsilon)/(sigma_sq1))/wzdeno -
      0.5 * sigx17_1)/sigx8_1 - sigx2 * sigx7_1 * sqrt(sigma_sq1)/sigx8_1^2) *
      ewu1 * ewz), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * ((sigx5_2 * dmusig2 *
      depsi2 + (S * (sigx15_2 - S * sigx5_2 * dmusig2 *
      sigx16_2 * (epsilon)) * (epsilon) - 0.5 * sigx17_2)/(sigma_sq2))/sigx8_2 -
      sigx2 * sigx7_2 * sqrt(sigma_sq2)/sigx8_2^2) * prC *
      ewu2), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((S * (sigx15_1 + S * sigx9_1 * dmusig1 * sigx16_1 *
    ewu1 * (epsilon)/ssq1^2) * (epsilon)/(sigma_sq1) - sigx9_1 *
    dmusig1 * depsi1 * ewu1/ssq1^2)/wzdeno - 0.5 * sigx17_1)/sigx8_1 -
    sigx2 * sigx18_1 * sqrt(sigma_sq1)/sigx8_1^2) * ewv1 *
    ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((S * (sigx15_2 +
      S * s9sq2 * dmusig2 * sigx16_2 * ewu2 * (epsilon)/ssq2^2) *
      (epsilon) - 0.5 * sigx17_2)/(sigma_sq2) - s9sq2 *
      dmusig2 * depsi2 * ewu2/ssq2^2)/sigx8_2 - sigx2 *
      sigx10_2 * sqrt(sigma_sq2)/sigx8_2^2) * prC * ewv2),
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (s12sq1 * sigx1_1/(sigma_sq1)) -
      (sigx2 * sigx12/(2 * sigx3_2 + 2 * sigx3_1) + 2 *
        (prC * sigx1_2/(wzdeno * (sigma_sq2) * sqrt(sigma_sq2))))) *
      ewz/(2 * sigx3_2 + 2 * sigx3_1), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((ewu1 * (S * (0.5 * sigx24_1 - (0.5 * (S^2 * sigx5_1 *
    depsi1 * (epsilon)^2/(sigma_sq1)^2) - (((0.5 * (ewu1/(sigma_sq1)) +
    1 - 0.5 * (0.5 * prU1 + ewu1/(sigma_sq1))) * prU1 * ewv1/sigmastar1 +
    (2 - 2 * (usq1^2 * ewu1 * (sigma_sq1)/ssq1^2)) * sigmastar1)/ssq1^2 +
    S^2 * sigx5_1^2 * ewu1 * (epsilon)^2/ssq1) * depsi1) *
    dmusig1) * (epsilon)/wzdeno - 0.5 * sigx25_1) + S * sigx6_1 *
    (epsilon)/wzdeno - 0.5 * sigx11_1)/sigx8_1 - sigx28_1 *
    ewu1 * sigx7_1/sigx8_1^2) * ewu1 * ewz), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (prC * ewu1 * ewu2 * ewz * sigx7_1 *
      sigx7_2 * sqrt(sigma_sq2)/sigx27_1)), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU1 *
      ewv1) - S^2 * sigx9_1 * sigx5_1 * ewu1 * (epsilon)^2)/(sigma_sq1) +
      0.5 * ((ewu1/(sigma_sq1) - 1) * ewv1/(sigma_sq1) +
        1 - 0.5 * (prU1 * prV1))) * depsi1/sigmastar1 +
      0.5 * (S^2 * sigx9_1 * depsi1 * (epsilon)^2/(sigma_sq1)^2)) *
      ewu1 + sigx9_1 * (1 - 2 * (usq1 * ewu1 * (sigma_sq1) *
      sigmastar1/ssq1^2)) * depsi1) * dmusig1/ssq1^2 +
      0.5 * sigx24_1) * (epsilon)/wzdeno - 0.5 * sigx25_1)/sigx8_1 -
      sigx28_1 * sigx18_1/sigx8_1^2) * ewu1 * ewv1 * ewz),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (prC * ewu1 * ewv2 *
      ewz * sigx10_2 * sigx7_1 * sqrt(sigma_sq2)/sigx27_1)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * (sigx21_1 - S * s12sq1 * sigx5_1 * dmusig1 * depsi1 *
      (epsilon)) - 2 * (sigx12 * ewz * sigx7_1/sigx8_1)) *
    ewu1 * ewz/(2 * sigx3_2 + 2 * sigx3_1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((ewu2 * (S * (0.5 *
      sigx24_2 - (0.5 * (S^2 * sigx5_2 * depsi2 * (epsilon)^2/(sigma_sq2)^2) -
      (((0.5 * (ewu2/(sigma_sq2)) + 1 - 0.5 * (0.5 * prU2 +
        ewu2/(sigma_sq2))) * prU2 * ewv2/sigmastar2 +
        (2 - 2 * (usq2^2 * ewu2 * (sigma_sq2)/ssq2^2)) *
          sigmastar2)/ssq2^2 + S^2 * sigx5_2^2 * ewu2 *
        (epsilon)^2/ssq2) * depsi2) * dmusig2) * (epsilon) -
      0.5 * sigx25_2) + S * sigx6_2 * (epsilon) - 0.5 *
      s7sq2)/sigx8_2 - sigx28_2 * ewu2 * sigx7_2/sigx8_2^2) *
      prC * ewu2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (prC * ewu2 * ewv1 *
      ewz * sigx18_1 * sigx7_2 * sqrt(sigma_sq1)/(sigx8_1^2 *
      sqrt(sigma_sq2)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((S * (((((0.5 * (prU2 * ewv2) - S^2 * s9sq2 * sigx5_2 *
    ewu2 * (epsilon)^2)/(sigma_sq2) + 0.5 * ((ewu2/(sigma_sq2) -
    1) * ewv2/(sigma_sq2) + 1 - 0.5 * (prU2 * prV2))) * depsi2/sigmastar2 +
    0.5 * (S^2 * s9sq2 * depsi2 * (epsilon)^2/(sigma_sq2)^2)) *
    ewu2 + s9sq2 * (1 - 2 * (usq2 * ewu2 * (sigma_sq2) *
    sigmastar2/ssq2^2)) * depsi2) * dmusig2/ssq2^2 + 0.5 *
    sigx24_2) * (epsilon) - 0.5 * sigx25_2)/sigx8_2 - sigx28_2 *
    sigx10_2/sigx8_2^2) * prC * ewu2 * ewv2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (prC * (2 * (sigx12 * sigx7_2/(2 * sigx3_2 +
      2 * sigx3_1)) + 2 * (S * sigx6_2 * (epsilon)/wzdeno -
      0.5 * sigx27_2)) * ewu2 * ewz/sigx8_2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) -
      0.5 * (0.5 * prV1 + ewv1/(sigma_sq1))) * prV1 + S^2 *
      sigx9_1^2 * ewu1 * ewv1 * (epsilon)^2/(ssq1^2 * (sigma_sq1))) *
      depsi1 * ewu1/sigmastar1 + ((0.5 * (S^2 * depsi1 *
      (epsilon)^2/(sigma_sq1)^2) - 2 * (sigx9_1 * depsi1 *
      (sigma_sq1) * sigmastar1/ssq1^2)) * ewv1 + depsi1) *
      sigx9_1) * dmusig1 * ewu1/ssq1^2 + S * (0.5 * (ewv1 *
      (S * (sigx9_1 * dmusig1 * ewu1/ssq1^2 + 0.5 * sigx22_1) *
        (epsilon) - 2 * (pmusig1/(sigma_sq1)))) + 0.5 *
      pmusig1) * depsi1 * (epsilon)/(sigma_sq1)^2) * (epsilon)/wzdeno -
      (0.5 * (depsi1 * pmusig1) + 0.5 * (ewv1 * (S * sigx10_1 *
        (epsilon) - wzdeno^2 * depsi1 * pmusig1/wzxsq1^2))) *
        wzdeno/wzxsq1^2)/sigx8_1 - (0.5 * sigx26_1 +
      2 * (ewz * sigx18_1)) * ewv1 * sigx18_1/sigx8_1^2) *
      ewv1 * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (prC * ewv1 * ewv2 * ewz * sigx18_1 *
      sigx10_2 * sqrt(sigma_sq2)/sigx27_1)), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx21_1 + S * sigx9_1 *
      s12sq1 * dmusig1 * depsi1 * ewu1 * (epsilon)/ssq1^2) -
      2 * (sigx12 * ewz * sigx18_1/sigx8_1)) * ewv1 * ewz/(2 *
      sigx3_2 + 2 * sigx3_1), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv2/(sigma_sq2)) -
      0.5 * (0.5 * prV2 + ewv2/(sigma_sq2))) * prV2 + S^2 *
      s9sq2^2 * ewu2 * ewv2 * (epsilon)^2/(ssq2^2 * (sigma_sq2))) *
      depsi2 * ewu2/sigmastar2 + ((0.5 * (S^2 * depsi2 *
      (epsilon)^2/(sigma_sq2)^2) - 2 * (s9sq2 * depsi2 *
      (sigma_sq2) * sigmastar2/ssq2^2)) * ewv2 + depsi2) *
      s9sq2) * dmusig2 * ewu2/ssq2^2 + S * (0.5 * (ewv2 *
      (S * (s9sq2 * dmusig2 * ewu2/ssq2^2 + 0.5 * sigx22_2) *
        (epsilon) - 2 * (pmusig2/(sigma_sq2)))) + 0.5 *
      pmusig2) * depsi2 * (epsilon)/(sigma_sq2)^2) * (epsilon) -
      (0.5 * (depsi2 * pmusig2) + 0.5 * (ewv2 * (S * sigx9_2 *
        (epsilon) - depsi2 * pmusig2/(sigma_sq2))))/(sigma_sq2))/sigx8_2 -
      (0.5 * sigx26_2 + 2 * (prC * sigx10_2)) * ewv2 *
        sigx10_2/sigx8_2^2) * prC * ewv2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (prC * (2 * (sigx12 * sigx10_2/(2 *
      sigx3_2 + 2 * sigx3_1)) + 2 * (S * sigx9_2 * (epsilon)/wzdeno -
      0.5 * sigx27_2)) * ewv2 * ewz/sigx8_2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    ((2 * (prC * (1/(wzdeno^2 * sqrt(sigma_sq2)) + sqrt(sigma_sq2)/(wzdeno *
      sqrt(sigma_sq2))^2) * depsi2 * pmusig2) - (sigx12^2/(2 *
      sigx3_2 + 2 * sigx3_1) + 2 * ((2 - 2 * (wzdeno *
      (sigma_sq1) * ewz/wzxsq1^2)) * depsi1 * pmusig1 *
      sqrt(sigma_sq1)/wzxsq1^2))) * ewz + 2 * (s12sq1 *
      depsi1 * pmusig1) - 2 * sigx11_2) * ewz/(2 * sigx3_2 +
    2 * sigx3_1), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf halfnormal-normal distribution
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
cnsfhalfnormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfhalfnormlike(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(ccnsfhalfnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfhalfnormlike, grad = cgradcnsfhalfnormlike,
      hess = chesscnsfhalfnormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfhalfnormlike(mleObj$par,
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
      mleObj$hessian <- chesscnsfhalfnormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfhalfnormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfhalfnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfhalfnormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

# different sigma_u
mcesfhalfnormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfhalfnormlike(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(cmcesfhalfnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfhalfnormlike, grad = cgradmcesfhalfnormlike,
      hess = chessmcesfhalfnormlike, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfhalfnormlike(mleObj$par,
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
      mleObj$hessian <- chessmcesfhalfnormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfhalfnormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfhalfnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfhalfnormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd

# same sigma_u
ccnsfhalfnormeff <- function(object, level) {
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
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
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
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) *
      pnorm(mustar1/sigmastar1 + sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) *
      pnorm(mustar2/sigmastar2 + sigmastar2)/pnorm(mustar2/sigmastar2)
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
cmcesfhalfnormeff <- function(object, level) {
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
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
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
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) *
      pnorm(mustar1/sigmastar1 + sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) *
      pnorm(mustar2/sigmastar2 + sigmastar2)/pnorm(mustar2/sigmastar2)
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
#' marginal impact on efficiencies for cnsf halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
# same sigma_u
ccnsfmarghalfnorm_Eu <- function(object) {
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
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmarghalfnorm_Vu <- function(object) {
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
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# different sigma_u
cmcesfmarghalfnorm_Eu <- function(object) {
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
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1/2) * dnorm(0), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2/2) * dnorm(0), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmarghalfnorm_Vu <- function(object) {
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
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}