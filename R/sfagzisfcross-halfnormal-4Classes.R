################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Generalized Zero Inefficiency Stochastic Frontier Analysis            #
# Number of Classes: 4L                                                        #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for gzisf 4 classes halfnormal-normal distribution
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
cGZISFhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(S * epsilon3/sqrt(exp(Wu3) + exp(Wv3))) *
    pnorm(mustar3/sigmastar3)
  Pi4 <- 1/exp(Wv4/2) * dnorm(S * epsilon4/exp(Wv4/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 + Probc4 * Pi4 <= 0, return(-.Machine$double.xmax),
    return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 + Probc4 *
      Pi4)))
}

# starting value for the log-likelihood ----------
#' starting values for gzisf 4 classes halfnormal-normal distribution
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
csGZISFfhalfnorm4C <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter, initAlg, printInfo,
  tol) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA + halfnormal - normal distributions...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = printInfo,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:(nXvar + 1)], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar +
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:nXvar], Esti[nXvar +
    1], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.98 * Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar -
    1), 0.98 * Esti[1:nXvar], Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar -
    1), rep(0, 3 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), names(Esti)[1:nXvar], paste0("Zv_", colnames(vHvar)), paste0("Cl1_",
    colnames(Zvar)), paste0("Cl2_", colnames(Zvar)), paste0("Cl3_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for gzisf 4 classes halfnormal-normal distribution
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
cgradGZISFhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewu3 <- exp(Wu3)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv3 <- exp(Wv3)
  ewv4_h <- exp(Wv4/2)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  ewz3 <- exp(Wz3)
  wzdeno <- (1 + ewz1 + ewz2 + ewz3)
  prC <- (1 - (ewz1 + ewz2 + ewz3)/wzdeno)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigma_sq3 <- ewu3 + ewv3
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  sigmastar3 <- sqrt(ewu3 * ewv3/(sigma_sq3))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  ssq3 <- ((sigma_sq3) * sigmastar3)
  musig1 <- (S * ewu1 * (epsilon1)/ssq1)
  musig2 <- (S * ewu2 * (epsilon2)/ssq2)
  musig3 <- (S * ewu3 * (epsilon3)/ssq3)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  dmusig3 <- dnorm(-musig3, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  pmusig3 <- pnorm(-musig3)
  epsi1 <- S * (epsilon1)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon2)/sqrt(sigma_sq2)
  epsi3 <- S * (epsilon3)/sqrt(sigma_sq3)
  depsi1 <- dnorm(epsi1)
  depsi2 <- dnorm(epsi2)
  depsi3 <- dnorm(epsi3)
  epsi4 <- S * (epsilon4)/ewv4_h
  depsi4 <- dnorm(epsi4)
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  prV3 <- (1 - ewv3/(sigma_sq3))
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  prU3 <- (1 - ewu3/(sigma_sq3))
  puv1 <- (prU1 * ewv1/sigmastar1)
  puv2 <- (prU2 * ewv2/sigmastar2)
  puv3 <- (prU3 * ewv3/sigmastar3)
  pvu1 <- (prV1 * ewu1/sigmastar1)
  pvu2 <- (prV2 * ewu2/sigmastar2)
  pvu3 <- (prV3 * ewu3/sigmastar3)
  dpesq1 <- (depsi1 * pmusig1/(sigma_sq1))
  dpesq2 <- (depsi2 * pmusig2/(sigma_sq2))
  dpesq3 <- (depsi3 * pmusig3/(sigma_sq3))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 * pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 * pmusig2 * (epsilon2))
  sigx1_3 <- (dmusig3 * depsi3 * ewu3/sigmastar3 + S * depsi3 * pmusig3 * (epsilon3))
  sigx2_1 <- (depsi1 * ewz1 * pmusig1/sqrt(sigma_sq1))
  sigx2_2 <- (depsi2 * ewz2 * pmusig2/sqrt(sigma_sq2))
  sigx2_3 <- (depsi3 * ewz3 * pmusig3/sqrt(sigma_sq3))
  sigx3 <- (prC * depsi4/ewv4_h + (2 * sigx2_1 + 2 * sigx2_2 + 2 * sigx2_3)/wzdeno)
  sigx4_1 <- (sigx3 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sigx4_2 <- (sigx3 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2))
  sigx4_3 <- (sigx3 * wzdeno * (sigma_sq3) * sqrt(sigma_sq3))
  sigx5_1 <- (S * depsi1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx5_2 <- (S * depsi2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  sigx5_3 <- (S * depsi3 * pmusig3 * (epsilon3)/(sigma_sq3)^2)
  sigx6_1 <- (1/ssq1 - (0.5 * puv1 + sigmastar1) * ewu1/ssq1^2)
  sigx6_2 <- (1/ssq2 - (0.5 * puv2 + sigmastar2) * ewu2/ssq2^2)
  sigx6_3 <- (1/ssq3 - (0.5 * puv3 + sigmastar3) * ewu3/ssq3^2)
  sigx7_1 <- (0.5 * sigx5_1 - sigx6_1 * dmusig1 * depsi1)
  sigx7_2 <- (0.5 * sigx5_2 - sigx6_2 * dmusig2 * depsi2)
  sigx7_3 <- (0.5 * sigx5_3 - sigx6_3 * dmusig3 * depsi3)
  sigx8_1 <- (S * sigx7_1 * (epsilon1) - 0.5 * dpesq1)
  sigx8_2 <- (S * sigx7_2 * (epsilon2) - 0.5 * dpesq2)
  sigx8_3 <- (S * sigx7_3 * (epsilon3) - 0.5 * dpesq3)
  sigx9_1 <- (sigx3 * wzdeno * sqrt(sigma_sq1))
  sigx9_2 <- (sigx3 * wzdeno * sqrt(sigma_sq2))
  sigx9_3 <- (sigx3 * wzdeno * sqrt(sigma_sq3))
  sigx10_1 <- ((0.5 * pvu1 + sigmastar1) * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 *
    sigx5_1)
  sigx10_2 <- ((0.5 * pvu2 + sigmastar2) * dmusig2 * depsi2 * ewu2/ssq2^2 + 0.5 *
    sigx5_2)
  sigx10_3 <- ((0.5 * pvu3 + sigmastar3) * dmusig3 * depsi3 * ewu3/ssq3^2 + 0.5 *
    sigx5_3)
  sigx11_1 <- (S * sigx10_1 * (epsilon1) - 0.5 * dpesq1)
  sigx11_2 <- (S * sigx10_2 * (epsilon2) - 0.5 * dpesq2)
  sigx11_3 <- (S * sigx10_3 * (epsilon3) - 0.5 * dpesq3)
  sigx12_1 <- (2 * (depsi1 * pmusig1/sqrt(sigma_sq1)) - sigx3)
  sigx12_2 <- (2 * (depsi2 * pmusig2/sqrt(sigma_sq2)) - sigx3)
  sigx12_3 <- (2 * (depsi3 * pmusig3/sqrt(sigma_sq3)) - sigx3)
  sigx13 <- (0.5 * (S^2 * depsi4 * (epsilon4)^2/ewv4_h^2) - 0.5 * depsi4)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = 2 * (S * sigx1_1 * ewz1/sigx4_1),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (ewu1 * ewz1 * sigx8_1/sigx9_1),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (ewv1 * ewz1 * sigx11_1/sigx9_1),
    FUN = "*"), sweep(Xvar, MARGIN = 1, STATS = 2 * (S * sigx1_2 * ewz2/sigx4_2),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (ewu2 * ewz2 * sigx8_2/sigx9_2),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (ewv2 * ewz2 * sigx11_2/sigx9_2),
    FUN = "*"), sweep(Xvar, MARGIN = 1, STATS = 2 * (S * sigx1_3 * ewz3/sigx4_3),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (ewu3 * ewz3 * sigx8_3/sigx9_3),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (ewv3 * ewz3 * sigx11_3/sigx9_3),
    FUN = "*"), sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * depsi4 * (epsilon4)/(sigx3 *
    ewv4_h^3), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx13 * prC/(sigx3 *
    ewv4_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx12_1 * ewz1/(sigx3 *
    wzdeno), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx12_2 * ewz2/(sigx3 *
    wzdeno), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx12_3 * ewz3/(sigx3 *
    wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for gzisf 4 classes halfnormal-normal distribution
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
chessGZISFhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewu3 <- exp(Wu3)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv3 <- exp(Wv3)
  ewv4_h <- exp(Wv4/2)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  ewz3 <- exp(Wz3)
  wzdeno <- (1 + ewz1 + ewz2 + ewz3)
  prC <- (1 - (ewz1 + ewz2 + ewz3)/wzdeno)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigma_sq3 <- ewu3 + ewv3
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  sigmastar3 <- sqrt(ewu3 * ewv3/(sigma_sq3))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  ssq3 <- ((sigma_sq3) * sigmastar3)
  musig1 <- (S * ewu1 * (epsilon1)/ssq1)
  musig2 <- (S * ewu2 * (epsilon2)/ssq2)
  musig3 <- (S * ewu3 * (epsilon3)/ssq3)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  dmusig3 <- dnorm(-musig3, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  pmusig3 <- pnorm(-musig3)
  epsi1 <- S * (epsilon1)/sqrt(sigma_sq1)
  epsi2 <- S * (epsilon2)/sqrt(sigma_sq2)
  epsi3 <- S * (epsilon3)/sqrt(sigma_sq3)
  depsi1 <- dnorm(epsi1)
  depsi2 <- dnorm(epsi2)
  depsi3 <- dnorm(epsi3)
  epsi4 <- S * (epsilon4)/ewv4_h
  depsi4 <- dnorm(epsi4)
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  prV3 <- (1 - ewv3/(sigma_sq3))
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  prU3 <- (1 - ewu3/(sigma_sq3))
  puv1 <- (prU1 * ewv1/sigmastar1)
  puv2 <- (prU2 * ewv2/sigmastar2)
  puv3 <- (prU3 * ewv3/sigmastar3)
  pvu1 <- (prV1 * ewu1/sigmastar1)
  pvu2 <- (prV2 * ewu2/sigmastar2)
  pvu3 <- (prV3 * ewu3/sigmastar3)
  dpesq1 <- (depsi1 * pmusig1/(sigma_sq1))
  dpesq2 <- (depsi2 * pmusig2/(sigma_sq2))
  dpesq3 <- (depsi3 * pmusig3/(sigma_sq3))
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 * pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + S * depsi2 * pmusig2 * (epsilon2))
  sigx1_3 <- (dmusig3 * depsi3 * ewu3/sigmastar3 + S * depsi3 * pmusig3 * (epsilon3))
  sigx2_1 <- (depsi1 * ewz1 * pmusig1/sqrt(sigma_sq1))
  sigx2_2 <- (depsi2 * ewz2 * pmusig2/sqrt(sigma_sq2))
  sigx2_3 <- (depsi3 * ewz3 * pmusig3/sqrt(sigma_sq3))
  sigx3 <- (prC * depsi4/ewv4_h + (2 * sigx2_1 + 2 * sigx2_2 + 2 * sigx2_3)/wzdeno)
  sigx4_1 <- (sigx3 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sigx4_2 <- (sigx3 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2))
  sigx4_3 <- (sigx3 * wzdeno * (sigma_sq3) * sqrt(sigma_sq3))
  sigx5_1 <- (S * depsi1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx5_2 <- (S * depsi2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  sigx5_3 <- (S * depsi3 * pmusig3 * (epsilon3)/(sigma_sq3)^2)
  sigx6_1 <- (1/ssq1 - (0.5 * puv1 + sigmastar1) * ewu1/ssq1^2)
  sigx6_2 <- (1/ssq2 - (0.5 * puv2 + sigmastar2) * ewu2/ssq2^2)
  sigx6_3 <- (1/ssq3 - (0.5 * puv3 + sigmastar3) * ewu3/ssq3^2)
  sigx7_1 <- (0.5 * sigx5_1 - sigx6_1 * dmusig1 * depsi1)
  sigx7_2 <- (0.5 * sigx5_2 - sigx6_2 * dmusig2 * depsi2)
  sigx7_3 <- (0.5 * sigx5_3 - sigx6_3 * dmusig3 * depsi3)
  sigx8_1 <- (S * sigx7_1 * (epsilon1) - 0.5 * dpesq1)
  sigx8_2 <- (S * sigx7_2 * (epsilon2) - 0.5 * dpesq2)
  sigx8_3 <- (S * sigx7_3 * (epsilon3) - 0.5 * dpesq3)
  sigx9_1 <- (sigx3 * wzdeno * sqrt(sigma_sq1))
  sigx9_2 <- (sigx3 * wzdeno * sqrt(sigma_sq2))
  sigx9_3 <- (sigx3 * wzdeno * sqrt(sigma_sq3))
  sigx10_1 <- ((0.5 * pvu1 + sigmastar1) * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 *
    sigx5_1)
  sigx10_2 <- ((0.5 * pvu2 + sigmastar2) * dmusig2 * depsi2 * ewu2/ssq2^2 + 0.5 *
    sigx5_2)
  sigx10_3 <- ((0.5 * pvu3 + sigmastar3) * dmusig3 * depsi3 * ewu3/ssq3^2 + 0.5 *
    sigx5_3)
  sigx11_1 <- (S * sigx10_1 * (epsilon1) - 0.5 * dpesq1)
  sigx11_2 <- (S * sigx10_2 * (epsilon2) - 0.5 * dpesq2)
  sigx11_3 <- (S * sigx10_3 * (epsilon3) - 0.5 * dpesq3)
  sigx12_1 <- (2 * (depsi1 * pmusig1/sqrt(sigma_sq1)) - sigx3)
  sigx12_2 <- (2 * (depsi2 * pmusig2/sqrt(sigma_sq2)) - sigx3)
  sigx12_3 <- (2 * (depsi3 * pmusig3/sqrt(sigma_sq3)) - sigx3)
  sigx13 <- (0.5 * (S^2 * depsi4 * (epsilon4)^2/ewv4_h^2) - 0.5 * depsi4)
  sigx14_1 <- (S * (dmusig1 * ewu1/sigmastar1 + S * pmusig1 * (epsilon1)) * (epsilon1)/(sigma_sq1) -
    pmusig1)
  sigx14_2 <- (S * (dmusig2 * ewu2/sigmastar2 + S * pmusig2 * (epsilon2)) * (epsilon2)/(sigma_sq2) -
    pmusig2)
  sigx14_3 <- (S * (dmusig3 * ewu3/sigmastar3 + S * pmusig3 * (epsilon3)) * (epsilon3)/(sigma_sq3) -
    pmusig3)
  sigx15_1 <- (depsi1 * ewu1/ewv1 + depsi1)
  sigx15_2 <- (depsi2 * ewu2/ewv2 + depsi2)
  sigx15_3 <- (depsi3 * ewu3/ewv3 + depsi3)
  sigx16_1 <- (S * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx16_2 <- (S * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  sigx16_3 <- (S * pmusig3 * (epsilon3)/(sigma_sq3)^2)
  sigx17_1 <- (S * (0.5 * sigx16_1 - sigx6_1 * dmusig1) * (epsilon1) - 2 * (pmusig1/(sigma_sq1)))
  sigx17_2 <- (S * (0.5 * sigx16_2 - sigx6_2 * dmusig2) * (epsilon2) - 2 * (pmusig2/(sigma_sq2)))
  sigx17_3 <- (S * (0.5 * sigx16_3 - sigx6_3 * dmusig3) * (epsilon3) - 2 * (pmusig3/(sigma_sq3)))
  sigx18_1 <- ((2 - 2 * (ewz1/wzdeno))/(sigx3 * wzdeno) - 2 * (sigx12_1 * ewz1/(sigx3 *
    wzdeno)^2))
  sigx18_2 <- ((2 - 2 * (ewz2/wzdeno))/(sigx3 * wzdeno) - 2 * (sigx12_2 * ewz2/(sigx3 *
    wzdeno)^2))
  sigx18_3 <- (2 * (sigx12_3/(sigx3 * wzdeno)^2) + 2/(sigx3 * wzdeno^2))
  sigx19_1 <- ((sigx3 * ewv4_h)^2 * wzdeno * sqrt(sigma_sq1))
  sigx19_2 <- ((sigx3 * ewv4_h)^2 * wzdeno * sqrt(sigma_sq2))
  sigx19_3 <- ((sigx3 * ewv4_h)^2 * wzdeno * sqrt(sigma_sq3))
  sigx20 <- (sigx9_2^2 * (sigma_sq1) * sqrt(sigma_sq1))
  sigx21_1 <- ((sigma_sq1) * sqrt(sigma_sq1))
  sigx21_2 <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx21_3 <- (sigma_sq3) * sqrt(sigma_sq3)
  sigx22_1 <- (sigx3 * wzdeno/sqrt(sigma_sq1))
  sigx22_2 <- (sigx3 * wzdeno/sqrt(sigma_sq2))
  sigx22_3 <- (sigx3 * wzdeno/sqrt(sigma_sq3))
  sigx23_1 <- (0.5 * sigx22_1 + 2 * (ewz1 * sigx8_1))
  sigx23_2 <- (0.5 * sigx22_2 + 2 * (ewz2 * sigx8_2))
  sigx23_3 <- (0.5 * sigx22_3 + 2 * (ewz3 * sigx8_3))
  sigx24_1 <- (S * depsi1 * sigx17_1 * (epsilon1)/(sigma_sq1)^2)
  sigx24_2 <- (S * depsi2 * sigx17_2 * (epsilon2)/(sigma_sq2)^2)
  sigx24_3 <- (S * depsi3 * sigx17_3 * (epsilon3)/(sigma_sq3)^2)
  sigx25_1 <- ((sigx3 * wzdeno)^2 * sqrt(sigma_sq1))
  sigx25_2 <- ((sigx3 * wzdeno)^2 * sqrt(sigma_sq2))
  sigx25_3 <- ((sigx3 * ewv4_h^3)^2 * wzdeno * sqrt(sigma_sq3))
  sigx26_1 <- ((sigx3 * ewv4_h^3)^2 * wzdeno * sqrt(sigma_sq1))
  sigx26_2 <- ((sigx3 * ewv4_h^3)^2 * wzdeno * sqrt(sigma_sq2))
  sigx27_1 <- (2 * (sigx12_1/(sigx3 * wzdeno)^2) + 2/(sigx3 * wzdeno^2))
  sigx27_2 <- (2 * (sigx12_2/(sigx3 * wzdeno)^2) + 2/(sigx3 * wzdeno^2))
  sigx28_2 <- (sigx9_3^2 * (sigma_sq2) * sqrt(sigma_sq2))
  hessll <- matrix(nrow = 4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar, ncol = 4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S^2 * ((depsi1 * sigx14_1 + S * dmusig1 * sigx15_1 * ewu1 * (epsilon1)/ssq1)/sigx4_1 -
    2 * (sigx1_1^2 * ewz1/sigx4_1^2)) * ewz1), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S * ((sigx6_1 * dmusig1 * depsi1 + (S * ((0.5 * sigx14_1 -
      0.5 * pmusig1) * depsi1/(sigma_sq1) - S * sigx6_1 * dmusig1 * sigx15_1 *
      (epsilon1)) * (epsilon1) - 0.5 * (sigx1_1/(sigma_sq1)))/(sigma_sq1))/sigx9_1 -
      2 * (sigx1_1 * ewz1 * sigx8_1/(sigx9_1^2 * (sigma_sq1)))) * ewu1 * ewz1),
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((S * ((0.5 * sigx14_1 - 0.5 * pmusig1) *
      depsi1/(sigma_sq1) + S * (0.5 * pvu1 + sigmastar1) * dmusig1 * sigx15_1 *
      ewu1 * (epsilon1)/ssq1^2) * (epsilon1) - 0.5 * (sigx1_1/(sigma_sq1)))/(sigma_sq1) -
      (0.5 * pvu1 + sigmastar1) * dmusig1 * depsi1 * ewu1/ssq1^2)/sigx9_1 -
      2 * (sigx1_1 * ewz1 * sigx11_1/(sigx9_1^2 * (sigma_sq1)))) * ewv1 * ewz1),
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S^2 * sigx1_1 * sigx1_2 * (sigma_sq2) *
      ewz1 * ewz2 * sqrt(sigma_sq2)/(sigx4_2^2 * (sigma_sq1) * sqrt(sigma_sq1)))),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar * (4 * (S *
    sigx1_1 * ewu2 * ewz1 * ewz2 * sigx8_2 * sqrt(sigma_sq2)/sigx20)), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar * (4 *
    (S * sigx1_1 * ewv2 * ewz1 * ewz2 * sigx11_2 * sqrt(sigma_sq2)/sigx20)),
    FUN = "*"), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (4 * (S^2 * sigx1_1 * sigx1_3 * (sigma_sq3) * ewz1 * ewz3 * sqrt(sigma_sq3)/(sigx4_3^2 *
      (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (4 * (S * sigx1_1 * ewu3 * ewz1 * ewz3 * sigx8_3 * sqrt(sigma_sq3)/(sigx9_3^2 *
      (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (4 * (S * sigx1_1 * ewv3 * ewz1 * ewz3 * sigx11_3 * sqrt(sigma_sq3)/(sigx9_3^2 *
      (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (2 * (S^3 * prC * sigx1_1 * depsi4 * ewv4_h^3 * ewz1 * (epsilon4)/((sigx3 *
      ewv4_h^3)^2 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"),
    Xvar)
  hessll[1:nXvar, (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (2 * (S * sigx13 * prC * sigx1_1 * ewv4_h * ewz1/((sigx3 * ewv4_h)^2 * wzdeno *
      (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[1:nXvar, (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * sigx18_1 * sigx1_1 * ewz1/sigx21_1, FUN = "*"), Zvar)
  hessll[1:nXvar, (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = -wHvar * (S * sigx27_2 * sigx1_1 * ewz1 * ewz2/sigx21_1), FUN = "*"),
    Zvar)
  hessll[1:nXvar, (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (S * sigx18_3 * sigx1_1 * ewz1 * ewz3/sigx21_1),
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((ewu1 * (S * (0.5 * sigx24_1 - (0.5 * (S^2 *
      sigx6_1 * depsi1 * (epsilon1)^2/(sigma_sq1)^2) - (((0.5 * (ewu1/(sigma_sq1)) +
      1 - 0.5 * (0.5 * prU1 + ewu1/(sigma_sq1))) * prU1 * ewv1/sigmastar1 +
      (2 - 2 * ((0.5 * puv1 + sigmastar1)^2 * ewu1 * (sigma_sq1)/ssq1^2)) *
        sigmastar1)/ssq1^2 + S^2 * sigx6_1^2 * ewu1 * (epsilon1)^2/ssq1) *
      depsi1) * dmusig1) * (epsilon1) - 0.5 * ((S * sigx7_1 * (epsilon1) -
      depsi1 * pmusig1/(sigma_sq1))/(sigma_sq1))) + S * sigx7_1 * (epsilon1) -
      0.5 * dpesq1)/sigx9_1 - sigx23_1 * ewu1 * sigx8_1/sigx9_1^2) * ewu1 *
      ewz1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 2 * (((S *
    (((((0.5 * (prU1 * ewv1) - S^2 * (0.5 * pvu1 + sigmastar1) * sigx6_1 * ewu1 *
      (epsilon1)^2)/(sigma_sq1) + 0.5 * ((ewu1/(sigma_sq1) - 1) * ewv1/(sigma_sq1) +
      1 - 0.5 * (prU1 * prV1))) * depsi1/sigmastar1 + 0.5 * (S^2 * (0.5 * pvu1 +
      sigmastar1) * depsi1 * (epsilon1)^2/(sigma_sq1)^2)) * ewu1 + (0.5 * pvu1 +
      sigmastar1) * (1 - 2 * ((0.5 * puv1 + sigmastar1) * ewu1 * (sigma_sq1) *
      sigmastar1/ssq1^2)) * depsi1) * dmusig1/ssq1^2 + 0.5 * sigx24_1) * (epsilon1) -
    0.5 * ((S * sigx7_1 * (epsilon1) - depsi1 * pmusig1/(sigma_sq1))/(sigma_sq1)))/sigx9_1 -
    sigx23_1 * sigx11_1/sigx9_1^2) * ewu1 * ewv1 * ewz1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (4 * (S * sigx1_2 * ewu1 * (sigma_sq2) * ewz1 * ewz2 * sigx8_1 * sqrt(sigma_sq2)/(sigx4_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (4 * (ewu1 * ewu2 * ewz1 * ewz2 * sigx8_1 * sigx8_2 * sqrt(sigma_sq2)/(sigx9_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (ewu1 * ewv2 * ewz1 * ewz2 * sigx11_2 * sigx8_1 * sqrt(sigma_sq2)/(sigx9_2^2 *
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (S * sigx1_3 * ewu1 * (sigma_sq3) * ewz1 * ewz3 * sigx8_1 *
      sqrt(sigma_sq3)/(sigx4_3^2 * sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (ewu1 * ewu3 * ewz1 * ewz3 * sigx8_1 * sigx8_3 * sqrt(sigma_sq3)/(sigx9_3^2 *
      sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (ewu1 * ewv3 * ewz1 * ewz3 * sigx11_3 * sigx8_1 * sqrt(sigma_sq3)/(sigx9_3^2 *
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (2 * (S^2 * prC * depsi4 * ewu1 * ewv4_h^3 * ewz1 * sigx8_1 *
      (epsilon4)/sigx26_1)), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (2 * (sigx13 * prC * ewu1 * ewv4_h * ewz1 * sigx8_1/sigx19_1)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar +
    1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx18_1 * ewu1 * ewz1 * sigx8_1/sqrt(sigma_sq1),
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar +
    nZHvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx27_2 * ewu1 * ewz1 * ewz2 * sigx8_1/sqrt(sigma_sq1)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar +
    2 * nZHvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx18_3 * ewu1 * ewz1 * ewz3 * sigx8_1/sqrt(sigma_sq1)),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) - 0.5 * (0.5 * prV1 + ewv1/(sigma_sq1))) *
    prV1 + S^2 * (0.5 * pvu1 + sigmastar1)^2 * ewu1 * ewv1 * (epsilon1)^2/(ssq1^2 *
    (sigma_sq1))) * depsi1 * ewu1/sigmastar1 + ((0.5 * (S^2 * depsi1 * (epsilon1)^2/(sigma_sq1)^2) -
    2 * ((0.5 * pvu1 + sigmastar1) * depsi1 * (sigma_sq1) * sigmastar1/ssq1^2)) *
    ewv1 + depsi1) * (0.5 * pvu1 + sigmastar1)) * dmusig1 * ewu1/ssq1^2 + S *
    (0.5 * (ewv1 * (S * ((0.5 * pvu1 + sigmastar1) * dmusig1 * ewu1/ssq1^2 +
      0.5 * sigx16_1) * (epsilon1) - 2 * (pmusig1/(sigma_sq1)))) + 0.5 * pmusig1) *
    depsi1 * (epsilon1)/(sigma_sq1)^2) * (epsilon1) - (0.5 * (depsi1 * pmusig1) +
    0.5 * (ewv1 * (S * sigx10_1 * (epsilon1) - depsi1 * pmusig1/(sigma_sq1))))/(sigma_sq1))/sigx9_1 -
    (0.5 * sigx22_1 + 2 * (ewz1 * sigx11_1)) * ewv1 * sigx11_1/sigx9_1^2) * ewv1 *
    ewz1), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (4 * (S * sigx1_2 * (sigma_sq2) * ewv1 * ewz1 * ewz2 * sigx11_1 *
      sqrt(sigma_sq2)/(sigx4_2^2 * sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (ewu2 * ewv1 * ewz1 * ewz2 * sigx11_1 *
      sigx8_2 * sqrt(sigma_sq2)/(sigx9_2^2 * sqrt(sigma_sq1)))), FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (ewv1 * ewv2 * ewz1 * ewz2 * sigx11_1 *
      sigx11_2 * sqrt(sigma_sq2)/(sigx9_2^2 * sqrt(sigma_sq1)))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S * sigx1_3 * (sigma_sq3) * ewv1 * ewz1 *
      ewz3 * sigx11_1 * sqrt(sigma_sq3)/(sigx4_3^2 * sqrt(sigma_sq1)))), FUN = "*"),
    Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (ewu3 * ewv1 * ewz1 * ewz3 * sigx11_1 *
      sigx8_3 * sqrt(sigma_sq3)/(sigx9_3^2 * sqrt(sigma_sq1)))), FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (4 * (ewv1 * ewv3 * ewz1 * ewz3 * sigx11_1 *
      sigx11_3 * sqrt(sigma_sq3)/(sigx9_3^2 * sqrt(sigma_sq1)))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^2 * prC * depsi4 * ewv1 * ewv4_h^3 *
      ewz1 * sigx11_1 * (epsilon4)/sigx26_1)), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx13 * prC * ewv1 * ewv4_h * ewz1 *
      sigx11_1/sigx19_1)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx18_1 * ewv1 * ewz1 * sigx11_1/sqrt(sigma_sq1),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx27_2 * ewv1 * ewz1 * ewz2 * sigx11_1/sqrt(sigma_sq1)),
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 *
    nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (sigx18_3 *
    ewv1 * ewz1 * ewz3 * sigx11_1/sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S^2 * ((depsi2 * sigx14_2 + S * dmusig2 *
      sigx15_2 * ewu2 * (epsilon2)/ssq2)/sigx4_2 - 2 * (sigx1_2^2 * ewz2/sigx4_2^2)) *
      ewz2), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * ((sigx6_2 * dmusig2 * depsi2 + (S *
      ((0.5 * sigx14_2 - 0.5 * pmusig2) * depsi2/(sigma_sq2) - S * sigx6_2 *
        dmusig2 * sigx15_2 * (epsilon2)) * (epsilon2) - 0.5 * (sigx1_2/(sigma_sq2)))/(sigma_sq2))/sigx9_2 -
      2 * (sigx1_2 * ewz2 * sigx8_2/(sigx9_2^2 * (sigma_sq2)))) * ewu2 * ewz2),
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((S * ((0.5 * sigx14_2 - 0.5 * pmusig2) *
      depsi2/(sigma_sq2) + S * (0.5 * pvu2 + sigmastar2) * dmusig2 * sigx15_2 *
      ewu2 * (epsilon2)/ssq2^2) * (epsilon2) - 0.5 * (sigx1_2/(sigma_sq2)))/(sigma_sq2) -
      (0.5 * pvu2 + sigmastar2) * dmusig2 * depsi2 * ewu2/ssq2^2)/sigx9_2 -
      2 * (sigx1_2 * ewz2 * sigx11_2/(sigx9_2^2 * (sigma_sq2)))) * ewv2 * ewz2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S^2 * sigx1_2 * sigx1_3 * (sigma_sq3) *
      ewz2 * ewz3 * sqrt(sigma_sq3)/(sigx4_3^2 * (sigma_sq2) * sqrt(sigma_sq2)))),
    FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S * sigx1_2 * ewu3 * ewz2 * ewz3 * sigx8_3 *
      sqrt(sigma_sq3)/sigx28_2)), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (4 * (S * sigx1_2 * ewv3 * ewz2 * ewz3 * sigx11_3 *
      sqrt(sigma_sq3)/sigx28_2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^3 * prC * sigx1_2 * depsi4 * ewv4_h^3 *
      ewz2 * (epsilon4)/((sigx3 * ewv4_h^3)^2 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2)))),
    FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S * sigx13 * prC * sigx1_2 * ewv4_h *
      ewz2/((sigx3 * ewv4_h)^2 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2)))),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar * (S * sigx27_1 *
    sigx1_2 * ewz1 * ewz2/sigx21_2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * sigx18_2 * sigx1_2 * ewz2/sigx21_2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (S * sigx18_3 * sigx1_2 * ewz2 * ewz3/sigx21_2), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((ewu2 * (S * (0.5 * sigx24_2 - (0.5 * (S^2 *
      sigx6_2 * depsi2 * (epsilon2)^2/(sigma_sq2)^2) - (((0.5 * (ewu2/(sigma_sq2)) +
      1 - 0.5 * (0.5 * prU2 + ewu2/(sigma_sq2))) * prU2 * ewv2/sigmastar2 +
      (2 - 2 * ((0.5 * puv2 + sigmastar2)^2 * ewu2 * (sigma_sq2)/ssq2^2)) *
        sigmastar2)/ssq2^2 + S^2 * sigx6_2^2 * ewu2 * (epsilon2)^2/ssq2) *
      depsi2) * dmusig2) * (epsilon2) - 0.5 * ((S * sigx7_2 * (epsilon2) -
      depsi2 * pmusig2/(sigma_sq2))/(sigma_sq2))) + S * sigx7_2 * (epsilon2) -
      0.5 * dpesq2)/sigx9_2 - sigx23_2 * ewu2 * sigx8_2/sigx9_2^2) * ewu2 *
      ewz2), FUN = "*"), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU2 * ewv2) - S^2 * (0.5 *
      pvu2 + sigmastar2) * sigx6_2 * ewu2 * (epsilon2)^2)/(sigma_sq2) + 0.5 *
      ((ewu2/(sigma_sq2) - 1) * ewv2/(sigma_sq2) + 1 - 0.5 * (prU2 * prV2))) *
      depsi2/sigmastar2 + 0.5 * (S^2 * (0.5 * pvu2 + sigmastar2) * depsi2 *
      (epsilon2)^2/(sigma_sq2)^2)) * ewu2 + (0.5 * pvu2 + sigmastar2) * (1 -
      2 * ((0.5 * puv2 + sigmastar2) * ewu2 * (sigma_sq2) * sigmastar2/ssq2^2)) *
      depsi2) * dmusig2/ssq2^2 + 0.5 * sigx24_2) * (epsilon2) - 0.5 * ((S *
      sigx7_2 * (epsilon2) - depsi2 * pmusig2/(sigma_sq2))/(sigma_sq2)))/sigx9_2 -
      sigx23_2 * sigx11_2/sigx9_2^2) * ewu2 * ewv2 * ewz2), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (4 *
    (S * sigx1_3 * ewu2 * (sigma_sq3) * ewz2 * ewz3 * sigx8_2 * sqrt(sigma_sq3)/(sigx4_3^2 *
      sqrt(sigma_sq2)))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (4 *
    (ewu2 * ewu3 * ewz2 * ewz3 * sigx8_2 * sigx8_3 * sqrt(sigma_sq3)/(sigx9_3^2 *
      sqrt(sigma_sq2)))), FUN = "*"), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (4 *
    (ewu2 * ewv3 * ewz2 * ewz3 * sigx11_3 * sigx8_2 * sqrt(sigma_sq3)/(sigx9_3^2 *
      sqrt(sigma_sq2)))), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (2 *
    (S^2 * prC * depsi4 * ewu2 * ewv4_h^3 * ewz2 * sigx8_2 * (epsilon4)/sigx26_2)),
    FUN = "*"), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 4 *
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (2 *
    (sigx13 * prC * ewu2 * ewv4_h * ewz2 * sigx8_2/sigx19_2)), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 4 *
      nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx27_1 * ewu2 * ewz1 * ewz2 * sigx8_2/sqrt(sigma_sq2)), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 3 * nuZUvar +
      4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx18_2 * ewu2 * ewz2 * sigx8_2/sqrt(sigma_sq2), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 3 *
      nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx18_3 * ewu2 * ewz2 * ewz3 * sigx8_2/sqrt(sigma_sq2)),
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar * 2 * (((S *
    ((((0.5 * (ewv2/(sigma_sq2)) - 0.5 * (0.5 * prV2 + ewv2/(sigma_sq2))) * prV2 +
      S^2 * (0.5 * pvu2 + sigmastar2)^2 * ewu2 * ewv2 * (epsilon2)^2/(ssq2^2 *
        (sigma_sq2))) * depsi2 * ewu2/sigmastar2 + ((0.5 * (S^2 * depsi2 *
      (epsilon2)^2/(sigma_sq2)^2) - 2 * ((0.5 * pvu2 + sigmastar2) * depsi2 *
      (sigma_sq2) * sigmastar2/ssq2^2)) * ewv2 + depsi2) * (0.5 * pvu2 + sigmastar2)) *
      dmusig2 * ewu2/ssq2^2 + S * (0.5 * (ewv2 * (S * ((0.5 * pvu2 + sigmastar2) *
      dmusig2 * ewu2/ssq2^2 + 0.5 * sigx16_2) * (epsilon2) - 2 * (pmusig2/(sigma_sq2)))) +
      0.5 * pmusig2) * depsi2 * (epsilon2)/(sigma_sq2)^2) * (epsilon2) - (0.5 *
    (depsi2 * pmusig2) + 0.5 * (ewv2 * (S * sigx10_2 * (epsilon2) - depsi2 *
    pmusig2/(sigma_sq2))))/(sigma_sq2))/sigx9_2 - (0.5 * sigx22_2 + 2 * (ewz2 *
    sigx11_2)) * ewv2 * sigx11_2/sigx9_2^2) * ewv2 * ewz2), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (4 *
    (S * sigx1_3 * (sigma_sq3) * ewv2 * ewz2 * ewz3 * sigx11_2 * sqrt(sigma_sq3)/(sigx4_3^2 *
      sqrt(sigma_sq2)))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (4 *
    (ewu3 * ewv2 * ewz2 * ewz3 * sigx11_2 * sigx8_3 * sqrt(sigma_sq3)/(sigx9_3^2 *
      sqrt(sigma_sq2)))), FUN = "*"), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (4 *
    (ewv2 * ewv3 * ewz2 * ewz3 * sigx11_2 * sigx11_3 * sqrt(sigma_sq3)/(sigx9_3^2 *
      sqrt(sigma_sq2)))), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (2 *
    (S^2 * prC * depsi4 * ewv2 * ewv4_h^3 * ewz2 * sigx11_2 * (epsilon4)/sigx26_2)),
    FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (2 *
    (sigx13 * prC * ewv2 * ewv4_h * ewz2 * sigx11_2/sigx19_2)), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx27_1 * ewv2 * ewz1 * ewz2 * sigx11_2/sqrt(sigma_sq2)), FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * sigx18_2 * ewv2 * ewz2 * sigx11_2/sqrt(sigma_sq2), FUN = "*"),
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx18_3 * ewv2 * ewz2 * ewz3 * sigx11_2/sqrt(sigma_sq2)),
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S^2 * ((depsi3 * sigx14_3 + S * dmusig3 * sigx15_3 * ewu3 * (epsilon3)/ssq3)/sigx4_3 -
    2 * (sigx1_3^2 * ewz3/sigx4_3^2)) * ewz3), FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * ((sigx6_3 * dmusig3 * depsi3 + (S * ((0.5 * sigx14_3 - 0.5 * pmusig3) *
    depsi3/(sigma_sq3) - S * sigx6_3 * dmusig3 * sigx15_3 * (epsilon3)) * (epsilon3) -
    0.5 * (sigx1_3/(sigma_sq3)))/(sigma_sq3))/sigx9_3 - 2 * (sigx1_3 * ewz3 *
    sigx8_3/(sigx9_3^2 * (sigma_sq3)))) * ewu3 * ewz3), FUN = "*"), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((S * ((0.5 * sigx14_3 - 0.5 * pmusig3) * depsi3/(sigma_sq3) +
    S * (0.5 * pvu3 + sigmastar3) * dmusig3 * sigx15_3 * ewu3 * (epsilon3)/ssq3^2) *
    (epsilon3) - 0.5 * (sigx1_3/(sigma_sq3)))/(sigma_sq3) - (0.5 * pvu3 + sigmastar3) *
    dmusig3 * depsi3 * ewu3/ssq3^2)/sigx9_3 - 2 * (sigx1_3 * ewz3 * sigx11_3/(sigx9_3^2 *
    (sigma_sq3)))) * ewv3 * ewz3), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (2 * (S^3 * prC * sigx1_3 * depsi4 * ewv4_h^3 * ewz3 * (epsilon4)/((sigx3 *
      ewv4_h^3)^2 * wzdeno * sigx21_3))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (2 * (S * sigx13 * prC * sigx1_3 * ewv4_h * ewz3/((sigx3 * ewv4_h)^2 * wzdeno *
      sigx21_3))), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (S * sigx27_1 * sigx1_3 * ewz1 * ewz3/(sigx21_3)), FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = -wHvar * (S * sigx27_2 * sigx1_3 * ewz2 * ewz3/(sigx21_3)), FUN = "*"),
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 - 2 * (ewz3/wzdeno))/(sigx3 * wzdeno) -
      2 * (sigx12_3 * ewz3/(sigx3 * wzdeno)^2)) * sigx1_3 * ewz3/(sigx21_3),
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((ewu3 * (S * (0.5 * sigx24_3 - (0.5 * (S^2 * sigx6_3 * depsi3 * (epsilon3)^2/(sigma_sq3)^2) -
    (((0.5 * (ewu3/(sigma_sq3)) + 1 - 0.5 * (0.5 * prU3 + ewu3/(sigma_sq3))) *
      prU3 * ewv3/sigmastar3 + (2 - 2 * ((0.5 * puv3 + sigmastar3)^2 * ewu3 *
      (sigma_sq3)/ssq3^2)) * sigmastar3)/ssq3^2 + S^2 * sigx6_3^2 * ewu3 *
      (epsilon3)^2/ssq3) * depsi3) * dmusig3) * (epsilon3) - 0.5 * ((S * sigx7_3 *
    (epsilon3) - depsi3 * pmusig3/(sigma_sq3))/(sigma_sq3))) + S * sigx7_3 *
    (epsilon3) - 0.5 * dpesq3)/sigx9_3 - sigx23_3 * ewu3 * sigx8_3/sigx9_3^2) *
    ewu3 * ewz3), FUN = "*"), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((S * (((((0.5 * (prU3 * ewv3) - S^2 * (0.5 * pvu3 + sigmastar3) * sigx6_3 *
    ewu3 * (epsilon3)^2)/(sigma_sq3) + 0.5 * ((ewu3/(sigma_sq3) - 1) * ewv3/(sigma_sq3) +
    1 - 0.5 * (prU3 * prV3))) * depsi3/sigmastar3 + 0.5 * (S^2 * (0.5 * pvu3 +
    sigmastar3) * depsi3 * (epsilon3)^2/(sigma_sq3)^2)) * ewu3 + (0.5 * pvu3 +
    sigmastar3) * (1 - 2 * ((0.5 * puv3 + sigmastar3) * ewu3 * (sigma_sq3) *
    sigmastar3/ssq3^2)) * depsi3) * dmusig3/ssq3^2 + 0.5 * sigx24_3) * (epsilon3) -
    0.5 * ((S * sigx7_3 * (epsilon3) - depsi3 * pmusig3/(sigma_sq3))/(sigma_sq3)))/sigx9_3 -
    sigx23_3 * sigx11_3/sigx9_3^2) * ewu3 * ewv3 * ewz3), FUN = "*"), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (S^2 * prC * depsi4 * ewu3 * ewv4_h^3 * ewz3 * sigx8_3 * (epsilon4)/sigx25_3)),
    FUN = "*"), Xvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (sigx13 * prC * ewu3 * ewv4_h * ewz3 * sigx8_3/sigx19_3)), FUN = "*"),
    vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (sigx27_1 * ewu3 * ewz1 * ewz3 * sigx8_3/sqrt(sigma_sq3)), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx27_2 * ewu3 * ewz2 * ewz3 * sigx8_3/sqrt(sigma_sq3)),
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 - 2 * (ewz3/wzdeno))/(sigx3 * wzdeno) - 2 *
      (sigx12_3 * ewz3/(sigx3 * wzdeno)^2)) * ewu3 * ewz3 * sigx8_3/sqrt(sigma_sq3),
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    2 * (((S * ((((0.5 * (ewv3/(sigma_sq3)) - 0.5 * (0.5 * prV3 + ewv3/(sigma_sq3))) *
    prV3 + S^2 * (0.5 * pvu3 + sigmastar3)^2 * ewu3 * ewv3 * (epsilon3)^2/(ssq3^2 *
    (sigma_sq3))) * depsi3 * ewu3/sigmastar3 + ((0.5 * (S^2 * depsi3 * (epsilon3)^2/(sigma_sq3)^2) -
    2 * ((0.5 * pvu3 + sigmastar3) * depsi3 * (sigma_sq3) * sigmastar3/ssq3^2)) *
    ewv3 + depsi3) * (0.5 * pvu3 + sigmastar3)) * dmusig3 * ewu3/ssq3^2 + S *
    (0.5 * (ewv3 * (S * ((0.5 * pvu3 + sigmastar3) * dmusig3 * ewu3/ssq3^2 +
      0.5 * sigx16_3) * (epsilon3) - 2 * (pmusig3/(sigma_sq3)))) + 0.5 * pmusig3) *
    depsi3 * (epsilon3)/(sigma_sq3)^2) * (epsilon3) - (0.5 * (depsi3 * pmusig3) +
    0.5 * (ewv3 * (S * sigx10_3 * (epsilon3) - depsi3 * pmusig3/(sigma_sq3))))/(sigma_sq3))/sigx9_3 -
    (0.5 * sigx22_3 + 2 * (ewz3 * sigx11_3)) * ewv3 * sigx11_3/sigx9_3^2) * ewv3 *
    ewz3), FUN = "*"), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (S^2 * prC * depsi4 * ewv3 * ewv4_h^3 * ewz3 * sigx11_3 * (epsilon4)/sigx25_3)),
    FUN = "*"), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (sigx13 * prC * ewv3 * ewv4_h * ewz3 * sigx11_3/sigx19_3)), FUN = "*"),
    vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx27_1 * ewv3 * ewz1 * ewz3 * sigx11_3/sqrt(sigma_sq3)), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx27_2 * ewv3 * ewz2 * ewz3 * sigx11_3/sqrt(sigma_sq3)),
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 - 2 * (ewz3/wzdeno))/(sigx3 * wzdeno) - 2 *
      (sigx12_3 * ewz3/(sigx3 * wzdeno)^2)) * ewv3 * ewz3 * sigx11_3/sqrt(sigma_sq3),
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * ((S^2 * (epsilon4)^2/ewv4_h^2 - 1)/(sigx3 * ewv4_h^3) - S^2 * prC *
    depsi4 * (epsilon4)^2/(sigx3 * ewv4_h^3)^2) * prC * depsi4, FUN = "*"), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * ((0.5 * (S^2 * (epsilon4)^2/ewv4_h^2 - 2) - 0.5)/(sigx3 * ewv4_h) -
    sigx13 * prC/(sigx3 * ewv4_h)^2) * prC * depsi4 * (epsilon4)/ewv4_h^2, FUN = "*"),
    vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
    (S^2 * (wzdeno * sigx12_1/(sigx3 * wzdeno)^2 + 1/(sigx3 * wzdeno)) * prC *
      depsi4 * ewz1 * (epsilon4)/ewv4_h^3), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = -wHvar * (S^2 * (wzdeno * sigx12_2/(sigx3 * wzdeno)^2 + 1/(sigx3 *
      wzdeno)) * prC * depsi4 * ewz2 * (epsilon4)/ewv4_h^3), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (S^2 * (wzdeno * sigx12_3/(sigx3 * wzdeno)^2 +
      1/(sigx3 * wzdeno)) * prC * depsi4 * ewz3 * (epsilon4)/ewv4_h^3), FUN = "*"),
    Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon4)^2/ewv4_h^2) - 1) - 0.25) * depsi4 *
    (epsilon4)^2/(sigx3 * ewv4_h^3) - (sigx13 * prC + 0.5 * (sigx3 * ewv4_h)) *
    sigx13/(sigx3 * ewv4_h)^2), FUN = "*"), vHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((wzdeno * sigx12_1/(sigx3 * wzdeno)^2 + 1/(sigx3 * wzdeno)) * sigx13 * prC *
      ewz1/ewv4_h), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * ((wzdeno * sigx12_2/(sigx3 * wzdeno)^2 + 1/(sigx3 * wzdeno)) *
      sigx13 * prC * ewz2/ewv4_h), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((wzdeno * sigx12_3/(sigx3 * wzdeno)^2 + 1/(sigx3 *
      wzdeno)) * sigx13 * prC * ewz3/ewv4_h), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * ((1 - ewz1/wzdeno)/(sigx3 * wzdeno) - 2 * (depsi1 * ewz1 *
      pmusig1/sigx25_1)) * sigx12_1 * ewz1, FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar +
    1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx12_1/(sigx3 * wzdeno^2) + 2 * (sigx12_2 *
      depsi1 * pmusig1/sigx25_1)) * ewz1 * ewz2), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    4 * nvZVvar + nZHvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar +
    1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx12_1/(sigx3 * wzdeno^2) + 2 * (sigx12_3 *
      depsi1 * pmusig1/sigx25_1)) * ewz1 * ewz3), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + 2 * nZHvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar +
    nZHvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((1 - ewz2/wzdeno)/(sigx3 * wzdeno) - 2 * (depsi2 *
      ewz2 * pmusig2/sigx25_2)) * sigx12_2 * ewz2, FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 3 *
    nuZUvar + 4 * nvZVvar + 2 * nZHvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar +
    2 * nZHvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx12_2/(sigx3 * wzdeno^2) + 2 * (sigx12_3 *
      depsi2 * pmusig2/sigx25_2)) * ewz2 * ewz3), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar +
    3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar), (4 * nXvar + 3 * nuZUvar + 4 * nvZVvar +
    2 * nZHvar + 1):(4 * nXvar + 3 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((1 - ewz3/wzdeno)/(sigx3 * wzdeno) - 2 * (depsi3 *
      ewz3 * pmusig3/((sigx3 * wzdeno)^2 * sqrt(sigma_sq3)))) * sigx12_3 *
      ewz3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for gzisf 4 classes halfnormal-normal distribution
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
GZISF4ChnormAlgOpt <- function(start, randStart, sdStart, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csGZISFfhalfnorm4C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, tol = tol,
      printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cGZISFhalfnormlike4C(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("GZISF 4 Classes Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cGZISFhalfnormlike4C(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradGZISFhalfnormlike4C(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0,
    control = list(trace = if (printInfo) 1 else 0, maxeval = itermax, stepmax = stepmax,
      xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cGZISFhalfnormlike4C,
    grad = cgradGZISFhalfnormlike4C, hess = chessGZISFhalfnormlike4C, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cGZISFhalfnormlike4C(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradGZISFhalfnormlike4C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cGZISFhalfnormlike4C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradGZISFhalfnormlike4C(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessGZISFhalfnormlike4C(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cGZISFhalfnormlike4C(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradGZISFhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar)), hess = function(parm) -chessGZISFhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar,
        nZHvar = nZHvar), print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cGZISFhalfnormlike4C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradGZISFhalfnormlike4C(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessGZISFhalfnormlike4C(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradGZISFhalfnormlike4C(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- chessGZISFhalfnormlike4C(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessGZISFhalfnormlike4C(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cGZISFhalfnormlike4C(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradGZISFhalfnormlike4C(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, if (is.null(start)) initHalf = initHalf))
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for gzisfcross 4 classes halfnormal-normal distribution
#' @param object object of class gzisfcross
#' @param level level for confidence interval
#' @noRd
cGZISF4Chalfnormeff <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(2 *
    object$nXvar + object$nuZUvar + object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    object$nZHvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    2 * object$nZHvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    3 * object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 1/exp(Wv4/2) * dnorm(object$S * epsilon4/exp(Wv4/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4), 1, which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c == 2, Pcond_c2, ifelse(Group_c ==
    3, Pcond_c3, Pcond_c4)))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
  u_c4 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, ifelse(Group_c == 2, u_c2, ifelse(Group_c ==
    3, u_c3, u_c4)))
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  ineff_c3 <- ifelse(Group_c == 3, u_c3, NA)
  ineff_c4 <- ifelse(Group_c == 4, u_c4, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c <- exp(-u_c)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c3 <- exp(-mustar3 + 1/2 * sigmastar3^2) * pnorm(mustar3/sigmastar3 -
      sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_c4 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, ifelse(Group_c == 2, teBC_c2, ifelse(Group_c ==
      3, teBC_c3, teBC_c4)))
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    effBC_c3 <- ifelse(Group_c == 3, teBC_c3, NA)
    effBC_c4 <- ifelse(Group_c == 4, teBC_c4, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c3 <- exp(mustar3 + 1/2 * sigmastar3^2) * pnorm(mustar3/sigmastar3 +
      sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_reciprocal_c4 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, ifelse(Group_c ==
      2, teBC_reciprocal_c2, ifelse(Group_c == 3, teBC_reciprocal_c3, teBC_reciprocal_c4)))
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    ReffBC_c3 <- ifelse(Group_c == 3, teBC_reciprocal_c3, NA)
    ReffBC_c4 <- ifelse(Group_c == 4, teBC_reciprocal_c4, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3, u_c3 = u_c3, teBC_c3 = teBC_c3,
      teBC_reciprocal_c3 = teBC_reciprocal_c3, PosteriorProb_c4 = Pcond_c4,
      PriorProb_c4 = Probc4, u_c4 = u_c4, teBC_c4 = teBC_c4, teBC_reciprocal_c4 = teBC_reciprocal_c4,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, ineff_c4 = ineff_c4,
      effBC_c1 = effBC_c1, effBC_c2 = effBC_c2, effBC_c3 = effBC_c3, effBC_c4 = effBC_c4,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2, ReffBC_c3 = ReffBC_c3,
      ReffBC_c4 = ReffBC_c4)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
      u_c3 = u_c3, PosteriorProb_c4 = Pcond_c4, PriorProb_c4 = Probc4, u_c4 = u_c4,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, ineff_c4 = ineff_c4)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal effects for for gzisf 4 classes halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargGZISF4Chalfnorm_Eu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(2 *
    object$nXvar + object$nuZUvar + object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    object$nZHvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    2 * object$nZHvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    3 * object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 1/exp(Wv4/2) * dnorm(object$S * epsilon4/exp(Wv4/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4), 1, which.max)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2) *
    dnorm(0), ncol = 1))
  colnames(margEff_c1) <- paste0("Eu_", colnames(uHvar)[-1], "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2) *
    dnorm(0), ncol = 1))
  colnames(margEff_c2) <- paste0("Eu_", colnames(uHvar)[-1], "_c2")
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar], nrow = 1), matrix(exp(Wu3/2) *
    dnorm(0), ncol = 1))
  colnames(margEff_c3) <- paste0("Eu_", colnames(uHvar)[-1], "_c3")
  margEff_c4 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  colnames(margEff_c4) <- paste0("Eu_", colnames(uHvar)[-1], "_c4")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_along(1:ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c], ifelse(Group_c ==
      2, margEff_c2[, c], ifelse(Group_c == 3, margEff_c3[, c], margEff_c4[,
      c])))
  }
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], "_c")
  margEff <- margEff <- data.frame(margEff_c, margEff_c1, margEff_c2, margEff_c3,
    margEff_c4)
  return(margEff)
}

cmargGZISF4Chalfnorm_Vu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(2 *
    object$nXvar + object$nuZUvar + object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    object$nZHvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    2 * object$nZHvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 4 * object$nvZVvar +
    3 * object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 1/exp(Wv4/2) * dnorm(object$S * epsilon4/exp(Wv4/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4), 1, which.max)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  colnames(margEff_c1) <- paste0("Vu_", colnames(uHvar)[-1], "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  colnames(margEff_c2) <- paste0("Vu_", colnames(uHvar)[-1], "_c2")
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar], nrow = 1), matrix(exp(Wu3) *
    (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  colnames(margEff_c3) <- paste0("Vu_", colnames(uHvar)[-1], "_c3")
  margEff_c4 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  colnames(margEff_c4) <- paste0("Vu_", colnames(uHvar)[-1], "_c4")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_along(1:ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c], ifelse(Group_c ==
      2, margEff_c2[, c], ifelse(Group_c == 3, margEff_c3[, c], margEff_c4[,
      c])))
  }
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2, margEff_c3, margEff_c4)
  return(margEff)
}
