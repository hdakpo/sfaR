################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Generalized Zero Inefficiency Stochastic Frontier Analysis            #
# Number of Classes: 3L                                                        #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for gzisf 3 classes halfnormal-normal distribution
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
cGZISFhalfnormlike3C <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar)]
  theta1 <- parm[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar)]
  theta2 <- parm[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 1/exp(Wv3/2) * dnorm(S * epsilon3/exp(Wv3/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2))
  Probc3 <- 1 - Probc1 - Probc2
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 <= 0, return(NA),
    return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2 + Probc3 *
      Pi3)))
}

# starting value for the log-likelihood ----------
#' starting values for gzisf 3 classes halfnormal-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param nZHvar number of separating variables
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
csGZISFfhalfnorm3C <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + halfnormal - normal distributions...\n")
  initHalf <- maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradhalfnormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = printInfo, reltol = tol), nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]),
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar)
  Esti <- initHalf$estimate
  StartVal <- c(Esti[1:(nXvar + 1)], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0,
    nvZVvar - 1), 0.98 * Esti[1:nXvar], Esti[nXvar + 1],
    if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar + 2],
    if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:nXvar],
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1),
    rep(0, 2 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    names(Esti)[1:nXvar], paste0("Zv_", colnames(vHvar)),
    paste0("Cl1_", colnames(Zvar)), paste0("Cl2_", colnames(Zvar)))
  names(initHalf$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for gzisf 3 classes halfnormal-normal distribution
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
cgradGZISFhalfnormlike3C <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar)]
  theta1 <- parm[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar)]
  theta2 <- parm[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv3_h <- exp(Wv3/2)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  wzdeno <- (1 + ewz1 + ewz2)
  prC <- (1 - (ewz1 + ewz2)/wzdeno)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu1 * (epsilon1)/ssq1)
  musig2 <- (S * ewu2 * (epsilon2)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsi1 <- S * (epsilon1)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon2)/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1)
  depsi2 <- dnorm(epsi2)
  epsi3 <- S * (epsilon3)/ewv3_h
  depsi3 <- dnorm(epsi3)
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  puv1 <- (prU1 * ewv1/sigmastar1)
  puv2 <- (prU2 * ewv2/sigmastar2)
  pvu1 <- (prV1 * ewu1/sigmastar1)
  pvu2 <- (prV2 * ewu2/sigmastar2)
  dpesq1 <- (depsi1 * pmusig1/(sigma_sq1))
  dpesq2 <- (depsi2 * pmusig2/(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 *
    pmusig2 * (epsilon2))
  sigx2_1 <- (depsi1 * ewz1 * pmusig1/sqrt(sigma_sq1))
  sigx2_2 <- (depsi2 * ewz2 * pmusig2/sqrt(sigma_sq2))
  sigx3 <- (prC * depsi3/ewv3_h + (2 * sigx2_1 + 2 * sigx2_2)/wzdeno)
  sigx4_1 <- (sigx3 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sigx4_2 <- (sigx3 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2))
  sigx5_1 <- (S * depsi1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx5_2 <- (S * depsi2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  sigx6_1 <- (1/ssq1 - (0.5 * puv1 + sigmastar1) * ewu1/ssq1^2)
  sigx6_2 <- (1/ssq2 - (0.5 * puv2 + sigmastar2) * ewu2/ssq2^2)
  sigx7_1 <- (0.5 * sigx5_1 - sigx6_1 * dmusig1 * depsi1)
  sigx7_2 <- (0.5 * sigx5_2 - sigx6_2 * dmusig2 * depsi2)
  sigx8_1 <- (S * sigx7_1 * (epsilon1) - 0.5 * dpesq1)
  sigx8_2 <- (S * sigx7_2 * (epsilon2) - 0.5 * dpesq2)
  sigx9_1 <- (sigx3 * wzdeno * sqrt(sigma_sq1))
  sigx9_2 <- (sigx3 * wzdeno * sqrt(sigma_sq2))
  sigx10_1 <- ((0.5 * pvu1 + sigmastar1) * dmusig1 * depsi1 *
    ewu1/ssq1^2 + 0.5 * sigx5_1)
  sigx10_2 <- ((0.5 * pvu2 + sigmastar2) * dmusig2 * depsi2 *
    ewu2/ssq2^2 + 0.5 * sigx5_2)
  sigx11_1 <- (S * sigx10_1 * (epsilon1) - 0.5 * dpesq1)
  sigx11_2 <- (S * sigx10_2 * (epsilon2) - 0.5 * dpesq2)
  sigx12_1 <- (2 * (depsi1 * pmusig1/sqrt(sigma_sq1)) - sigx3)
  sigx12_2 <- (2 * (depsi2 * pmusig2/sqrt(sigma_sq2)) - sigx3)
  sigx13 <- (0.5 * (S^2 * depsi3 * (epsilon3)^2/ewv3_h^2) -
    0.5 * depsi3)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = 2 * (S *
    sigx1_1 * ewz1/sigx4_1), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = 2 * (ewu1 * ewz1 * sigx8_1/sigx9_1), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (ewv1 * ewz1 * sigx11_1/sigx9_1),
      FUN = "*"), sweep(Xvar, MARGIN = 1, STATS = 2 * (S *
      sigx1_2 * ewz2/sigx4_2), FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = 2 * (ewu2 * ewz2 * sigx8_2/sigx9_2),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
      (ewv2 * ewz2 * sigx11_2/sigx9_2), FUN = "*"), sweep(Xvar,
      MARGIN = 1, STATS = S^2 * prC * depsi3 * (epsilon3)/(sigx3 *
        ewv3_h^3), FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx13 * prC/(sigx3 * ewv3_h), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx12_1 * ewz1/(sigx3 *
      wzdeno), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx12_2 *
      ewz2/(sigx3 * wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for gzisf 3 classes halfnormal-normal distribution
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
chessGZISFhalfnormlike3C <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar)]
  theta1 <- parm[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar)]
  theta2 <- parm[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv3_h <- exp(Wv3/2)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  wzdeno <- (1 + ewz1 + ewz2)
  prC <- (1 - (ewz1 + ewz2)/wzdeno)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  musig1 <- (S * ewu1 * (epsilon1)/ssq1)
  musig2 <- (S * ewu2 * (epsilon2)/ssq2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsi1 <- S * (epsilon1)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon2)/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1)
  depsi2 <- dnorm(epsi2)
  epsi3 <- S * (epsilon3)/ewv3_h
  depsi3 <- dnorm(epsi3)
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  puv1 <- (prU1 * ewv1/sigmastar1)
  puv2 <- (prU2 * ewv2/sigmastar2)
  pvu1 <- (prV1 * ewu1/sigmastar1)
  pvu2 <- (prV2 * ewu2/sigmastar2)
  dpesq1 <- (depsi1 * pmusig1/(sigma_sq1))
  dpesq2 <- (depsi2 * pmusig2/(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 *
    pmusig2 * (epsilon2))
  sigx2_1 <- (depsi1 * ewz1 * pmusig1/sqrt(sigma_sq1))
  sigx2_2 <- (depsi2 * ewz2 * pmusig2/sqrt(sigma_sq2))
  sigx3 <- (prC * depsi3/ewv3_h + (2 * sigx2_1 + 2 * sigx2_2)/wzdeno)
  sigx4_1 <- (sigx3 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sigx4_2 <- (sigx3 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2))
  sigx5_1 <- (S * depsi1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx5_2 <- (S * depsi2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  sigx6_1 <- (1/ssq1 - (0.5 * puv1 + sigmastar1) * ewu1/ssq1^2)
  sigx6_2 <- (1/ssq2 - (0.5 * puv2 + sigmastar2) * ewu2/ssq2^2)
  sigx7_1 <- (0.5 * sigx5_1 - sigx6_1 * dmusig1 * depsi1)
  sigx7_2 <- (0.5 * sigx5_2 - sigx6_2 * dmusig2 * depsi2)
  sigx8_1 <- (S * sigx7_1 * (epsilon1) - 0.5 * dpesq1)
  sigx8_2 <- (S * sigx7_2 * (epsilon2) - 0.5 * dpesq2)
  sigx9_1 <- (sigx3 * wzdeno * sqrt(sigma_sq1))
  sigx9_2 <- (sigx3 * wzdeno * sqrt(sigma_sq2))
  sigx10_1 <- ((0.5 * pvu1 + sigmastar1) * dmusig1 * depsi1 *
    ewu1/ssq1^2 + 0.5 * sigx5_1)
  sigx10_2 <- ((0.5 * pvu2 + sigmastar2) * dmusig2 * depsi2 *
    ewu2/ssq2^2 + 0.5 * sigx5_2)
  sigx11_1 <- (S * sigx10_1 * (epsilon1) - 0.5 * dpesq1)
  sigx11_2 <- (S * sigx10_2 * (epsilon2) - 0.5 * dpesq2)
  sigx12_1 <- (2 * (depsi1 * pmusig1/sqrt(sigma_sq1)) - sigx3)
  sigx12_2 <- (2 * (depsi2 * pmusig2/sqrt(sigma_sq2)) - sigx3)
  sigx13 <- (0.5 * (S^2 * depsi3 * (epsilon3)^2/ewv3_h^2) -
    0.5 * depsi3)
  sigx14_1 <- (S * (dmusig1 * ewu1/sigmastar1 + S * pmusig1 *
    (epsilon1)) * (epsilon1)/(sigma_sq1) - pmusig1)
  sigx14_2 <- (S * (dmusig2 * ewu2/sigmastar2 + S * pmusig2 *
    (epsilon2)) * (epsilon2)/(sigma_sq2) - pmusig2)
  sigx15_1 <- (depsi1 * ewu1/ewv1 + depsi1)
  sigx15_2 <- (depsi2 * ewu2/ewv2 + depsi2)
  sigx16_1 <- (S * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx16_2 <- (S * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  sigx17_1 <- (S * (0.5 * sigx16_1 - sigx6_1 * dmusig1) * (epsilon1) -
    2 * (pmusig1/(sigma_sq1)))
  sigx17_2 <- (S * (0.5 * sigx16_2 - sigx6_2 * dmusig2) * (epsilon2) -
    2 * (pmusig2/(sigma_sq2)))
  sigx18_1 <- ((2 - 2 * (ewz1/wzdeno))/(sigx3 * wzdeno) - 2 *
    (sigx12_1 * ewz1/(sigx3 * wzdeno)^2))
  sigx18_2 <- (2 * (sigx12_2/(sigx3 * wzdeno)^2) + 2/(sigx3 *
    wzdeno^2))
  sigx19_1 <- ((sigx3 * ewv3_h)^2 * wzdeno * sqrt(sigma_sq1))
  sigx19_2 <- ((sigx3 * ewv3_h)^2 * wzdeno * sqrt(sigma_sq2))
  sigx20 <- (sigx9_2^2 * (sigma_sq1) * sqrt(sigma_sq1))
  sigx21_1 <- ((sigma_sq1) * sqrt(sigma_sq1))
  sigx21_2 <- (sigma_sq2) * sqrt(sigma_sq2)
  sigx22_1 <- (sigx3 * wzdeno/sqrt(sigma_sq1))
  sigx22_2 <- (sigx3 * wzdeno/sqrt(sigma_sq2))
  sigx23_1 <- (0.5 * sigx22_1 + 2 * (ewz1 * sigx8_1))
  sigx23_2 <- (0.5 * sigx22_2 + 2 * (ewz2 * sigx8_2))
  sigx24_1 <- (S * depsi1 * sigx17_1 * (epsilon1)/(sigma_sq1)^2)
  sigx24_2 <- (S * depsi2 * sigx17_2 * (epsilon2)/(sigma_sq2)^2)
  sigx25_1 <- ((sigx3 * wzdeno)^2 * sqrt(sigma_sq1))
  sigx25_2 <- ((sigx3 * ewv3_h^3)^2 * wzdeno * sqrt(sigma_sq2))
  sigx26_1 <- ((sigx3 * ewv3_h^3)^2 * wzdeno * sqrt(sigma_sq1))
  sigx27_1 <- (2 * (sigx12_1/(sigx3 * wzdeno)^2) + 2/(sigx3 *
    wzdeno^2))
  sigx27_2 <- ((2 - 2 * (ewz2/wzdeno))/(sigx3 * wzdeno) - 2 *
    (sigx12_2 * ewz2/(sigx3 * wzdeno)^2))
  hessll <- matrix(nrow = 3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    2 * nZHvar, ncol = 3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    2 * nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S^2 * ((depsi1 * sigx14_1 + S *
      dmusig1 * sigx15_1 * ewu1 * (epsilon1)/ssq1)/sigx4_1 -
      2 * (sigx1_1^2 * ewz1/sigx4_1^2)) * ewz1), FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * ((sigx6_1 * dmusig1 *
      depsi1 + (S * ((0.5 * sigx14_1 - 0.5 * pmusig1) *
      depsi1/(sigma_sq1) - S * sigx6_1 * dmusig1 * sigx15_1 *
      (epsilon1)) * (epsilon1) - 0.5 * (sigx1_1/(sigma_sq1)))/(sigma_sq1))/sigx9_1 -
      2 * (sigx1_1 * ewz1 * sigx8_1/(sigx9_1^2 * (sigma_sq1)))) *
      ewu1 * ewz1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((S * ((0.5 * sigx14_1 - 0.5 * pmusig1) * depsi1/(sigma_sq1) +
    S * (0.5 * pvu1 + sigmastar1) * dmusig1 * sigx15_1 *
      ewu1 * (epsilon1)/ssq1^2) * (epsilon1) - 0.5 * (sigx1_1/(sigma_sq1)))/(sigma_sq1) -
    (0.5 * pvu1 + sigmastar1) * dmusig1 * depsi1 * ewu1/ssq1^2)/sigx9_1 -
    2 * (sigx1_1 * ewz1 * sigx11_1/(sigx9_1^2 * (sigma_sq1)))) *
    ewv1 * ewz1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = -wHvar * (4 * (S^2 * sigx1_1 * sigx1_2 * (sigma_sq2) *
      ewz1 * ewz2 * sqrt(sigma_sq2)/(sigx4_2^2 * (sigma_sq1) *
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S * sigx1_1 * ewu2 *
      ewz1 * ewz2 * sigx8_2 * sqrt(sigma_sq2)/sigx20)),
    FUN = "*"), uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S * sigx1_1 * ewv2 *
      ewz1 * ewz2 * sigx11_2 * sqrt(sigma_sq2)/sigx20)),
    FUN = "*"), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^3 * prC * sigx1_1 *
      depsi3 * ewv3_h^3 * ewz1 * (epsilon3)/((sigx3 * ewv3_h^3)^2 *
      wzdeno * (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"),
    Xvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S * sigx13 * prC *
      sigx1_1 * ewv3_h * ewz1/((sigx3 * ewv3_h)^2 * wzdeno *
      (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * sigx18_1 * sigx1_1 *
      ewz1/sigx21_1, FUN = "*"), Zvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    nZHvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (S * sigx18_2 * sigx1_1 * ewz1 * ewz2/sigx21_1), FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((ewu1 * (S * (0.5 * sigx24_1 - (0.5 * (S^2 * sigx6_1 *
    depsi1 * (epsilon1)^2/(sigma_sq1)^2) - (((0.5 * (ewu1/(sigma_sq1)) +
    1 - 0.5 * (0.5 * prU1 + ewu1/(sigma_sq1))) * prU1 * ewv1/sigmastar1 +
    (2 - 2 * ((0.5 * puv1 + sigmastar1)^2 * ewu1 * (sigma_sq1)/ssq1^2)) *
      sigmastar1)/ssq1^2 + S^2 * sigx6_1^2 * ewu1 * (epsilon1)^2/ssq1) *
    depsi1) * dmusig1) * (epsilon1) - 0.5 * ((S * sigx7_1 *
    (epsilon1) - depsi1 * pmusig1/(sigma_sq1))/(sigma_sq1))) +
    S * sigx7_1 * (epsilon1) - 0.5 * dpesq1)/sigx9_1 - sigx23_1 *
    ewu1 * sigx8_1/sigx9_1^2) * ewu1 * ewz1), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU1 *
      ewv1) - S^2 * (0.5 * pvu1 + sigmastar1) * sigx6_1 *
      ewu1 * (epsilon1)^2)/(sigma_sq1) + 0.5 * ((ewu1/(sigma_sq1) -
      1) * ewv1/(sigma_sq1) + 1 - 0.5 * (prU1 * prV1))) *
      depsi1/sigmastar1 + 0.5 * (S^2 * (0.5 * pvu1 + sigmastar1) *
      depsi1 * (epsilon1)^2/(sigma_sq1)^2)) * ewu1 + (0.5 *
      pvu1 + sigmastar1) * (1 - 2 * ((0.5 * puv1 + sigmastar1) *
      ewu1 * (sigma_sq1) * sigmastar1/ssq1^2)) * depsi1) *
      dmusig1/ssq1^2 + 0.5 * sigx24_1) * (epsilon1) - 0.5 *
      ((S * sigx7_1 * (epsilon1) - depsi1 * pmusig1/(sigma_sq1))/(sigma_sq1)))/sigx9_1 -
      sigx23_1 * sigx11_1/sigx9_1^2) * ewu1 * ewv1 * ewz1),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S * sigx1_2 * ewu1 *
      (sigma_sq2) * ewz1 * ewz2 * sigx8_1 * sqrt(sigma_sq2)/(sigx4_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (ewu1 * ewu2 * ewz1 *
      ewz2 * sigx8_1 * sigx8_2 * sqrt(sigma_sq2)/(sigx9_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (ewu1 * ewv2 * ewz1 *
      ewz2 * sigx11_2 * sigx8_1 * sqrt(sigma_sq2)/(sigx9_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^2 * prC * depsi3 *
      ewu1 * ewv3_h^3 * ewz1 * sigx8_1 * (epsilon3)/sigx26_1)),
    FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx13 * prC * ewu1 *
      ewv3_h * ewz1 * sigx8_1/sigx19_1)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx18_1 * ewu1 * ewz1 * sigx8_1/sqrt(sigma_sq1), FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + nZHvar + 1):(3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx18_2 * ewu1 * ewz1 *
      ewz2 * sigx8_1/sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) -
      0.5 * (0.5 * prV1 + ewv1/(sigma_sq1))) * prV1 + S^2 *
      (0.5 * pvu1 + sigmastar1)^2 * ewu1 * ewv1 * (epsilon1)^2/(ssq1^2 *
      (sigma_sq1))) * depsi1 * ewu1/sigmastar1 + ((0.5 *
      (S^2 * depsi1 * (epsilon1)^2/(sigma_sq1)^2) - 2 *
      ((0.5 * pvu1 + sigmastar1) * depsi1 * (sigma_sq1) *
        sigmastar1/ssq1^2)) * ewv1 + depsi1) * (0.5 *
      pvu1 + sigmastar1)) * dmusig1 * ewu1/ssq1^2 + S *
      (0.5 * (ewv1 * (S * ((0.5 * pvu1 + sigmastar1) *
        dmusig1 * ewu1/ssq1^2 + 0.5 * sigx16_1) * (epsilon1) -
        2 * (pmusig1/(sigma_sq1)))) + 0.5 * pmusig1) *
      depsi1 * (epsilon1)/(sigma_sq1)^2) * (epsilon1) -
      (0.5 * (depsi1 * pmusig1) + 0.5 * (ewv1 * (S * sigx10_1 *
        (epsilon1) - depsi1 * pmusig1/(sigma_sq1))))/(sigma_sq1))/sigx9_1 -
      (0.5 * sigx22_1 + 2 * (ewz1 * sigx11_1)) * ewv1 *
        sigx11_1/sigx9_1^2) * ewv1 * ewz1), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (4 * (S * sigx1_2 * (sigma_sq2) * ewv1 * ewz1 * ewz2 *
      sigx11_1 * sqrt(sigma_sq2)/(sigx4_2^2 * sqrt(sigma_sq1)))),
    FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 *
      nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (ewu2 * ewv1 * ewz1 * ewz2 * sigx11_1 *
      sigx8_2 * sqrt(sigma_sq2)/(sigx9_2^2 * sqrt(sigma_sq1)))),
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar +
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (ewv1 * ewv2 * ewz1 *
      ewz2 * sigx11_1 * sigx11_2 * sqrt(sigma_sq2)/(sigx9_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^2 * prC * depsi3 *
      ewv1 * ewv3_h^3 * ewz1 * sigx11_1 * (epsilon3)/sigx26_1)),
    FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
      2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx13 * prC * ewv1 *
      ewv3_h * ewz1 * sigx11_1/sigx19_1)), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar +
      2 * nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx18_1 * ewv1 * ewz1 *
      sigx11_1/sqrt(sigma_sq1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 *
      nXvar + 2 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx18_2 * ewv1 * ewz1 *
      ewz2 * sigx11_1/sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S^2 * ((depsi2 * sigx14_2 + S *
      dmusig2 * sigx15_2 * ewu2 * (epsilon2)/ssq2)/sigx4_2 -
      2 * (sigx1_2^2 * ewz2/sigx4_2^2)) * ewz2), FUN = "*"),
    Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * ((sigx6_2 * dmusig2 * depsi2 +
      (S * ((0.5 * sigx14_2 - 0.5 * pmusig2) * depsi2/(sigma_sq2) -
        S * sigx6_2 * dmusig2 * sigx15_2 * (epsilon2)) *
        (epsilon2) - 0.5 * (sigx1_2/(sigma_sq2)))/(sigma_sq2))/sigx9_2 -
      2 * (sigx1_2 * ewz2 * sigx8_2/(sigx9_2^2 * (sigma_sq2)))) *
      ewu2 * ewz2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((S * ((0.5 * sigx14_2 -
      0.5 * pmusig2) * depsi2/(sigma_sq2) + S * (0.5 *
      pvu2 + sigmastar2) * dmusig2 * sigx15_2 * ewu2 *
      (epsilon2)/ssq2^2) * (epsilon2) - 0.5 * (sigx1_2/(sigma_sq2)))/(sigma_sq2) -
      (0.5 * pvu2 + sigmastar2) * dmusig2 * depsi2 * ewu2/ssq2^2)/sigx9_2 -
      2 * (sigx1_2 * ewz2 * sigx11_2/(sigx9_2^2 * (sigma_sq2)))) *
      ewv2 * ewz2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^3 * prC * sigx1_2 *
      depsi3 * ewv3_h^3 * ewz2 * (epsilon3)/((sigx3 * ewv3_h^3)^2 *
      wzdeno * sigx21_2))), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S * sigx13 * prC *
      sigx1_2 * ewv3_h * ewz2/((sigx3 * ewv3_h)^2 * wzdeno *
      sigx21_2))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 1):(3 *
    nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (S * sigx27_1 * sigx1_2 *
      ewz1 * ewz2/(sigx21_2)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx27_2 * sigx1_2 *
      ewz2/(sigx21_2)), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar), (2 * nXvar + nuZUvar + nvZVvar +
    1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((ewu2 * (S * (0.5 *
      sigx24_2 - (0.5 * (S^2 * sigx6_2 * depsi2 * (epsilon2)^2/(sigma_sq2)^2) -
      (((0.5 * (ewu2/(sigma_sq2)) + 1 - 0.5 * (0.5 * prU2 +
        ewu2/(sigma_sq2))) * prU2 * ewv2/sigmastar2 +
        (2 - 2 * ((0.5 * puv2 + sigmastar2)^2 * ewu2 *
          (sigma_sq2)/ssq2^2)) * sigmastar2)/ssq2^2 +
        S^2 * sigx6_2^2 * ewu2 * (epsilon2)^2/ssq2) *
        depsi2) * dmusig2) * (epsilon2) - 0.5 * ((S *
      sigx7_2 * (epsilon2) - depsi2 * pmusig2/(sigma_sq2))/(sigma_sq2))) +
      S * sigx7_2 * (epsilon2) - 0.5 * dpesq2)/sigx9_2 -
      (0.5 * sigx22_2 + 2 * (ewz2 * sigx8_2)) * ewu2 *
        sigx8_2/sigx9_2^2) * ewu2 * ewz2), FUN = "*"),
    uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar +
    1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU2 *
      ewv2) - S^2 * (0.5 * pvu2 + sigmastar2) * sigx6_2 *
      ewu2 * (epsilon2)^2)/(sigma_sq2) + 0.5 * ((ewu2/(sigma_sq2) -
      1) * ewv2/(sigma_sq2) + 1 - 0.5 * (prU2 * prV2))) *
      depsi2/sigmastar2 + 0.5 * (S^2 * (0.5 * pvu2 + sigmastar2) *
      depsi2 * (epsilon2)^2/(sigma_sq2)^2)) * ewu2 + (0.5 *
      pvu2 + sigmastar2) * (1 - 2 * ((0.5 * puv2 + sigmastar2) *
      ewu2 * (sigma_sq2) * sigmastar2/ssq2^2)) * depsi2) *
      dmusig2/ssq2^2 + 0.5 * sigx24_2) * (epsilon2) - 0.5 *
      ((S * sigx7_2 * (epsilon2) - depsi2 * pmusig2/(sigma_sq2))/(sigma_sq2)))/sigx9_2 -
      (0.5 * sigx22_2 + 2 * (ewz2 * sigx8_2)) * sigx11_2/sigx9_2^2) *
      ewu2 * ewv2 * ewz2), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^2 * prC * depsi3 *
      ewu2 * ewv3_h^3 * ewz2 * sigx8_2 * (epsilon3)/sigx25_2)),
    FUN = "*"), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx13 * prC * ewu2 *
      ewv3_h * ewz2 * sigx8_2/sigx19_2)), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx27_1 * ewu2 * ewz1 *
      ewz2 * sigx8_2/sqrt(sigma_sq2)), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    nZHvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    2 * nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (sigx27_2 * ewu2 * ewz2 * sigx8_2/sqrt(sigma_sq2)), FUN = "*"),
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv2/(sigma_sq2)) -
      0.5 * (0.5 * prV2 + ewv2/(sigma_sq2))) * prV2 + S^2 *
      (0.5 * pvu2 + sigmastar2)^2 * ewu2 * ewv2 * (epsilon2)^2/(ssq2^2 *
      (sigma_sq2))) * depsi2 * ewu2/sigmastar2 + ((0.5 *
      (S^2 * depsi2 * (epsilon2)^2/(sigma_sq2)^2) - 2 *
      ((0.5 * pvu2 + sigmastar2) * depsi2 * (sigma_sq2) *
        sigmastar2/ssq2^2)) * ewv2 + depsi2) * (0.5 *
      pvu2 + sigmastar2)) * dmusig2 * ewu2/ssq2^2 + S *
      (0.5 * (ewv2 * (S * ((0.5 * pvu2 + sigmastar2) *
        dmusig2 * ewu2/ssq2^2 + 0.5 * sigx16_2) * (epsilon2) -
        2 * (pmusig2/(sigma_sq2)))) + 0.5 * pmusig2) *
      depsi2 * (epsilon2)/(sigma_sq2)^2) * (epsilon2) -
      (0.5 * (depsi2 * pmusig2) + 0.5 * (ewv2 * (S * sigx10_2 *
        (epsilon2) - depsi2 * pmusig2/(sigma_sq2))))/(sigma_sq2))/sigx9_2 -
      (0.5 * sigx22_2 + 2 * (ewz2 * sigx11_2)) * ewv2 *
        sigx11_2/sigx9_2^2) * ewv2 * ewz2), FUN = "*"),
    vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^2 * prC * depsi3 *
      ewv2 * ewv3_h^3 * ewz2 * sigx11_2 * (epsilon3)/sigx25_2)),
    FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx13 * prC * ewv2 *
      ewv3_h * ewz2 * sigx11_2/sigx19_2)), FUN = "*"),
    vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx27_1 * ewv2 * ewz1 * ewz2 * sigx11_2/sqrt(sigma_sq2)),
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + nZHvar + 1):(3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx27_2 * ewv2 * ewz2 *
      sigx11_2/sqrt(sigma_sq2), FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * ((S^2 * (epsilon3)^2/ewv3_h^2 -
      1)/(sigx3 * ewv3_h^3) - S^2 * prC * depsi3 * (epsilon3)^2/(sigx3 *
      ewv3_h^3)^2) * prC * depsi3, FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * ((0.5 * (S^2 * (epsilon3)^2/ewv3_h^2 -
      2) - 0.5)/(sigx3 * ewv3_h) - sigx13 * prC/(sigx3 *
      ewv3_h)^2) * prC * depsi3 * (epsilon3)/ewv3_h^2,
    FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (S^2 * (wzdeno * sigx12_1/(sigx3 * wzdeno)^2 + 1/(sigx3 *
      wzdeno)) * prC * depsi3 * ewz1 * (epsilon3)/ewv3_h^3),
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + nZHvar + 1):(3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = -wHvar * (S^2 * (wzdeno * sigx12_2/(sigx3 * wzdeno)^2 +
      1/(sigx3 * wzdeno)) * prC * depsi3 * ewz2 * (epsilon3)/ewv3_h^3),
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 3 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon3)^2/ewv3_h^2) - 1) - 0.25) * depsi3 *
      (epsilon3)^2/(sigx3 * ewv3_h^3) - (sigx13 * prC +
      0.5 * (sigx3 * ewv3_h)) * sigx13/(sigx3 * ewv3_h)^2),
    FUN = "*"), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 3 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar +
    nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((wzdeno * sigx12_1/(sigx3 * wzdeno)^2 + 1/(sigx3 * wzdeno)) *
      sigx13 * prC * ewz1/ewv3_h), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 3 * nvZVvar), (3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + nZHvar + 1):(3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((wzdeno * sigx12_2/(sigx3 *
      wzdeno)^2 + 1/(sigx3 * wzdeno)) * sigx13 * prC *
      ewz2/ewv3_h), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 3 * nvZVvar + nZHvar), (3 * nXvar + 2 *
    nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    3 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * ((1 - ewz1/wzdeno)/(sigx3 * wzdeno) -
      2 * (depsi1 * ewz1 * pmusig1/sigx25_1)) * sigx12_1 *
      ewz1, FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar +
    2 * nuZUvar + 3 * nvZVvar + nZHvar), (3 * nXvar + 2 *
    nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar + 2 *
    nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx12_1/(sigx3 * wzdeno^2) +
      2 * (sigx12_2 * depsi1 * pmusig1/sigx25_1)) * ewz1 *
      ewz2), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + 2 * nZHvar),
    (3 * nXvar + 2 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 *
      nXvar + 2 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((1 - ewz2/wzdeno)/(sigx3 *
      wzdeno) - 2 * (depsi2 * ewz2 * pmusig2/((sigx3 *
      wzdeno)^2 * sqrt(sigma_sq2)))) * sigx12_2 * ewz2,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for gzisf 3 classes halfnormal-normal distribution
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
GZISF3ChnormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csGZISFfhalfnorm3C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar,
      itermax = itermax, tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(cGZISFhalfnormlike3C(startVal, nXvar = nXvar,
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
  cat("GZISF 3 Classes Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cGZISFhalfnormlike3C(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradGZISFhalfnormlike3C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cGZISFhalfnormlike3C,
    grad = cgradGZISFhalfnormlike3C, hess = chessGZISFhalfnormlike3C,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cGZISFhalfnormlike3C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessGZISFhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradGZISFhalfnormlike3C(mleObj$par,
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
      mleObj$hessian <- chessGZISFhalfnormlike3C(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessGZISFhalfnormlike3C(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cGZISFhalfnormlike3C(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradGZISFhalfnormlike3C(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) initHalf = initHalf))
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for gzisf 3 classes halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
cGZISF3Chalfnormeff <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar)]
  theta1 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar + 1):(3 * object$nXvar +
    2 * object$nuZUvar + 3 * object$nvZVvar + 2 * object$nZHvar)]
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
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 1/exp(Wv3/2) * dnorm(object$S * epsilon3/exp(Wv3/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2))
  Probc3 <- 1 - Probc1 - Probc2
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3), 1,
    which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c ==
    2, Pcond_c2, Pcond_c3))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, ifelse(Group_c == 2, u_c2,
    u_c3))
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  ineff_c3 <- ifelse(Group_c == 3, u_c3, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c <- exp(-u_c)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c3 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, ifelse(Group_c ==
      2, teBC_c2, teBC_c3))
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    effBC_c3 <- ifelse(Group_c == 3, teBC_c3, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) *
      pnorm(mustar1/sigmastar1 + sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) *
      pnorm(mustar2/sigmastar2 + sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c3 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      ifelse(Group_c == 2, teBC_reciprocal_c2, teBC_reciprocal_c3))
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    ReffBC_c3 <- ifelse(Group_c == 3, teBC_reciprocal_c3,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, PosteriorProb_c3 = Pcond_c3,
      PriorProb_c3 = Probc3, u_c3 = u_c3, teBC_c3 = teBC_c3,
      teBC_reciprocal_c3 = teBC_reciprocal_c3, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, effBC_c3 = effBC_c3, ReffBC_c1 = ReffBC_c1,
      ReffBC_c2 = ReffBC_c2, ReffBC_c3 = ReffBC_c3)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
      u_c3 = u_c3, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2,
      ineff_c3 = ineff_c3)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal effects for for gzisf 3 classes halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargGZISF3Chalfnorm_Eu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar)]
  theta1 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar + 1):(3 * object$nXvar +
    2 * object$nuZUvar + 3 * object$nvZVvar + 2 * object$nZHvar)]
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
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 1/exp(Wv3/2) * dnorm(object$S * epsilon3/exp(Wv3/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2))
  Probc3 <- 1 - Probc1 - Probc2
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3), 1,
    which.max)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu1/2) * dnorm(0), ncol = 1))
  colnames(margEff_c1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu2/2) * dnorm(0), ncol = 1))
  colnames(margEff_c2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  margEff_c3 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  colnames(margEff_c3) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c3")
  margEff_c <- ifelse(Group_c == 1, margEff_c1, ifelse(Group_c ==
    2, margEff_c2, margEff_c3))
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  margEff <- bind_cols(margEff_c, margEff_c1, margEff_c2, margEff_c3)
  return(margEff)
}

cmargGZISF3Chalfnorm_Vu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar)]
  theta1 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar + 1):(3 * object$nXvar +
    2 * object$nuZUvar + 3 * object$nvZVvar + 2 * object$nZHvar)]
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
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 1/exp(Wv3/2) * dnorm(object$S * epsilon3/exp(Wv3/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2))
  Probc3 <- 1 - Probc1 - Probc2
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3), 1,
    which.max)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu1) * (1 - (dnorm(0)/pnorm(0))^2),
    ncol = 1))
  colnames(margEff_c1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu2) * (1 - (dnorm(0)/pnorm(0))^2),
    ncol = 1))
  colnames(margEff_c2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  margEff_c3 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  colnames(margEff_c3) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c3")
  margEff_c <- ifelse(Group_c == 1, margEff_c1, ifelse(Group_c ==
    2, margEff_c2, margEff_c3))
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  margEff <- bind_cols(margEff_c, margEff_c1, margEff_c2, margEff_c3)
  return(margEff)
}
