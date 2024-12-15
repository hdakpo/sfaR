################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Latent Class Stochastic Frontier Analysis                             #
# Number of Classes: 4L                                                        #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for lcm 4 classes halfnormal-normal distribution
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
cLCMhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 *
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
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
  mustar4 <- -exp(Wu4) * S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
 ll <- log(Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 + Probc4 * Pi4)
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for lcm 4 classes halfnormal-normal distribution
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
csLCMfhalfnorm4C <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA + halfnormal - normal distributions...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter,
        printLevel = printInfo, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)], Esti[nXvar +
    1], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar +
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)],
    Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1),
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1),
    0.98 * Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
      1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
      1) rep(0, nvZVvar - 1), rep(0, 3 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)),
    paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    paste0("Cl1_", colnames(Zvar)), paste0("Cl2_", colnames(Zvar)),
    paste0("Cl3_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for lcm 4 classes halfnormal-normal distribution
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
cgradLCMhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 *
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
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
  .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wu3)
  .e6 <- exp(Wv3)
  .e7 <- exp(Wu4)
  .e8 <- exp(Wv4)
  .e9 <- exp(Wz1)
  .e10 <- exp(Wz2)
  .e11 <- exp(Wz3)
  .e12 <- .e1 + .e2
  .e13 <- .e3 + .e4
  .e14 <- .e5 + .e6
  .e15 <- .e7 + .e8
  .e16 <- sqrt(.e1 * .e2/(.e12))
  .e17 <- sqrt(.e3 * .e4/(.e13))
  .e18 <- sqrt(.e5 * .e6/(.e14))
  .e19 <- sqrt(.e7 * .e8/(.e15))
  .e20 <- (S * .e1 * epsilon1/((.e12) * .e16))
  .e21 <- (S * .e3 * epsilon2/((.e13) * .e17))
  .e22 <- (S * .e5 * epsilon3/((.e14) * .e18))
  .e23 <- (S * .e7 * epsilon4/((.e15) * .e19))
  .e24 <- pnorm(-.e20)
  .e25 <- pnorm(-.e21)
  .e26 <- pnorm(-.e22)
  .e27 <- pnorm(-.e23)
  .e28 <- dnorm(-.e20)
  .e29 <- dnorm(S * epsilon1/sqrt(.e12))
  .e30 <- dnorm(S * epsilon2/sqrt(.e13))
  .e31 <- dnorm(S * epsilon3/sqrt(.e14))
  .e32 <- dnorm(S * epsilon4/sqrt(.e15))
  .e33 <- dnorm(-.e21)
  .e34 <- dnorm(-.e22)
  .e35 <- dnorm(-.e23)
  .e36 <- (.e29 * .e9 * .e24/sqrt(.e12))
  .e37 <- (.e30 * .e10 * .e25/sqrt(.e13))
  .e38 <- (.e31 * .e11 * .e26/sqrt(.e14))
  .e39 <- (.e28 * .e29 * .e1/.e16 + S * .e29 * .e24 * epsilon1)
  .e40 <- (1 - (.e9 + .e10 + .e11)/(1 + .e9 + .e10 + .e11))
  .e41 <- (.e40 * .e32 * .e27/sqrt(.e15))
  .e42 <- ((2 * .e36 + 2 * .e37 + 2 * .e38)/(1 + .e9 + .e10 + .e11) + 2 * .e41)
  .e43 <- (.e42 * (1 + .e9 + .e10 + .e11) * (.e12) * sqrt(.e12))
  .e44 <- (0.5 * ((1 - .e1/(.e12)) * .e2/.e16) + .e16)
  .e45 <- (1/((.e12) * .e16) - .e44 * .e1/((.e12) * .e16)^2)
  .e46 <- (S * .e29 * .e24 * epsilon1/(.e12)^2)
  .e47 <- (S * (0.5 * .e46 - .e45 * .e28 * .e29) * epsilon1 - 0.5 * (.e29 * .e24/(.e12)))
  .e48 <- (.e42 * (1 + .e9 + .e10 + .e11) * sqrt(.e12))
  .e49 <- (0.5 * ((1 - .e2/(.e12)) * .e1/.e16) + .e16)
  .e50 <- (.e49 * .e28 * .e29 * .e1/((.e12) * .e16)^2 + 0.5 * .e46)
  .e51 <- (S * .e50 * epsilon1 - 0.5 * (.e29 * .e24/(.e12)))
  .e52 <- (.e33 * .e30 * .e3/.e17 + S * .e30 * .e25 * epsilon2)
  .e53 <- (.e42 * (1 + .e9 + .e10 + .e11) * (.e13) * sqrt(.e13))
  .e54 <- (S * .e30 * .e25 * epsilon2/(.e13)^2)
  .e55 <- (0.5 * ((1 - .e3/(.e13)) * .e4/.e17) + .e17)
  .e56 <- (1/((.e13) * .e17) - .e55 * .e3/((.e13) * .e17)^2)
  .e57 <- (S * (0.5 * .e54 - .e56 * .e33 * .e30) * epsilon2 - 0.5 * (.e30 * .e25/(.e13)))
  .e58 <- (.e42 * (1 + .e9 + .e10 + .e11) * sqrt(.e13))
  .e59 <- (0.5 * ((1 - .e4/(.e13)) * .e3/.e17) + .e17)
  .e60 <- (.e59 * .e33 * .e30 * .e3/((.e13) * .e17)^2 + 0.5 * .e54)
  .e61 <- (S * .e60 * epsilon2 - 0.5 * (.e30 * .e25/(.e13)))
  .e62 <- (.e34 * .e31 * .e5/.e18 + S * .e31 * .e26 * epsilon3)
  .e63 <- (.e42 * (1 + .e9 + .e10 + .e11) * (.e14) * sqrt(.e14))
  .e64 <- (S * .e31 * .e26 * epsilon3/(.e14)^2)
  .e65 <- (0.5 * ((1 - .e5/(.e14)) * .e6/.e18) + .e18)
  .e66 <- (1/((.e14) * .e18) - .e65 * .e5/((.e14) * .e18)^2)
  .e67 <- (S * (0.5 * .e64 - .e66 * .e34 * .e31) * epsilon3 - 0.5 * (.e31 * .e26/(.e14)))
  .e68 <- (.e42 * (1 + .e9 + .e10 + .e11) * sqrt(.e14))
  .e69 <- (0.5 * ((1 - .e6/(.e14)) * .e5/.e18) + .e18)
  .e70 <- (.e69 * .e34 * .e31 * .e5/((.e14) * .e18)^2 + 0.5 * .e64)
  .e71 <- (S * .e70 * epsilon3 - 0.5 * (.e31 * .e26/(.e14)))
  .e72 <- (.e35 * .e32 * .e7/.e19 + S * .e32 * .e27 * epsilon4)
  .e73 <- (.e42 * (.e15) * sqrt(.e15))
  .e74 <- (S * .e32 * .e27 * epsilon4/(.e15)^2)
  .e75 <- (0.5 * ((1 - .e7/(.e15)) * .e8/.e19) + .e19)
  .e76 <- (1/((.e15) * .e19) - .e75 * .e7/((.e15) * .e19)^2)
  .e77 <- (S * (0.5 * .e74 - .e76 * .e35 * .e32) * epsilon4 - 0.5 * (.e32 * .e27/(.e15)))
  .e78 <- (0.5 * ((1 - .e8/(.e15)) * .e7/.e19) + .e19)
  .e79 <- (.e78 * .e35 * .e32 * .e7/((.e15) * .e19)^2 + 0.5 * .e74)
  .e80 <- (S * .e79 * epsilon4 - 0.5 * (.e32 * .e27/(.e15)))
  .e81 <- (2 * (.e29 * .e24/sqrt(.e12)) - .e42)
  .e82 <- (.e42 * (1 + .e9 + .e10 + .e11))
  .e83 <- (2 * (.e30 * .e25/sqrt(.e13)) - .e42)
  .e84 <- (2 * (.e31 * .e26/sqrt(.e14)) - .e42)
  gradll <- cbind(2 * (S * Xvar * .e39 * .e9/.e43), 2 * (uHvar * .e1 * .e9 * .e47/.e48),
    2 * (vHvar * .e2 * .e9 * .e51/.e48), 2 * (S * Xvar * .e52 * .e10/.e53), 2 *
      (uHvar * .e3 * .e10 * .e57/.e58), 2 * (vHvar * .e4 * .e10 * .e61/.e58),
    2 * (S * Xvar * .e62 * .e11/.e63), 2 * (uHvar * .e5 * .e11 * .e67/.e68),
    2 * (vHvar * .e6 * .e11 * .e71/.e68), 2 * (S * Xvar * .e40 * .e72/.e73),
    2 * (uHvar * .e40 * .e7 * .e77/(.e42 * sqrt(.e15))), 2 * (vHvar * .e40 *
      .e8 * .e80/(.e42 * sqrt(.e15))), Zvar * .e81 * .e9/.e82, Zvar * .e83 *
      .e10/.e82, Zvar * .e84 * .e11/.e82)
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for lcm 4 classes halfnormal-normal distribution
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
chessLCMhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 *
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
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
 .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wu3)
  .e6 <- exp(Wv3)
  .e7 <- exp(Wu4)
  .e8 <- exp(Wv4)
  .e9 <- exp(Wz1)
  .e10 <- exp(Wz2)
  .e11 <- exp(Wz3)
  .e12 <- .e1 + .e2
  .e13 <- .e3 + .e4
  .e14 <- .e5 + .e6
  .e15 <- .e7 + .e8
  .e16 <- sqrt(.e1 * .e2/(.e12))
  .e17 <- sqrt(.e3 * .e4/(.e13))
  .e18 <- sqrt(.e5 * .e6/(.e14))
  .e19 <- sqrt(.e7 * .e8/(.e15))
  .e20 <- (S * .e1 * epsilon1/((.e12) * .e16))
  .e21 <- (S * .e3 * epsilon2/((.e13) * .e17))
  .e22 <- (S * .e5 * epsilon3/((.e14) * .e18))
  .e23 <- (S * .e7 * epsilon4/((.e15) * .e19))
  .e24 <- pnorm(-.e20)
  .e25 <- pnorm(-.e21)
  .e26 <- pnorm(-.e22)
  .e27 <- pnorm(-.e23)
  .e28 <- dnorm(-.e20)
  .e29 <- dnorm(S * epsilon1/sqrt(.e12))
  .e30 <- dnorm(S * epsilon2/sqrt(.e13))
  .e31 <- dnorm(S * epsilon3/sqrt(.e14))
  .e32 <- dnorm(S * epsilon4/sqrt(.e15))
  .e33 <- dnorm(-.e21)
  .e34 <- dnorm(-.e22)
  .e35 <- dnorm(-.e23)
  .e36 <- (.e29 * .e9 * .e24/sqrt(.e12))
  .e37 <- (.e30 * .e10 * .e25/sqrt(.e13))
  .e38 <- (.e31 * .e11 * .e26/sqrt(.e14))
  .e39 <- (.e28 * .e29 * .e1/.e16 + S * .e29 * .e24 * epsilon1)
  .e40 <- (1 - (.e9 + .e10 + .e11)/(1 + .e9 + .e10 + .e11))
  .e41 <- (.e40 * .e32 * .e27/sqrt(.e15))
  .e42 <- ((2 * .e36 + 2 * .e37 + 2 * .e38)/(1 + .e9 + .e10 + .e11) + 2 * .e41)
  .e43 <- (.e42 * (1 + .e9 + .e10 + .e11) * (.e12) * sqrt(.e12))
  .e44 <- (0.5 * ((1 - .e1/(.e12)) * .e2/.e16) + .e16)
  .e45 <- (1/((.e12) * .e16) - .e44 * .e1/((.e12) * .e16)^2)
  .e46 <- (S * .e29 * .e24 * epsilon1/(.e12)^2)
  .e47 <- (S * (0.5 * .e46 - .e45 * .e28 * .e29) * epsilon1 - 0.5 * (.e29 * .e24/(.e12)))
  .e48 <- (.e42 * (1 + .e9 + .e10 + .e11) * sqrt(.e12))
  .e49 <- (0.5 * ((1 - .e2/(.e12)) * .e1/.e16) + .e16)
  .e50 <- (.e49 * .e28 * .e29 * .e1/((.e12) * .e16)^2 + 0.5 * .e46)
  .e51 <- (S * .e50 * epsilon1 - 0.5 * (.e29 * .e24/(.e12)))
  .e52 <- (.e33 * .e30 * .e3/.e17 + S * .e30 * .e25 * epsilon2)
  .e53 <- (.e42 * (1 + .e9 + .e10 + .e11) * (.e13) * sqrt(.e13))
  .e54 <- (S * .e30 * .e25 * epsilon2/(.e13)^2)
  .e55 <- (0.5 * ((1 - .e3/(.e13)) * .e4/.e17) + .e17)
  .e56 <- (1/((.e13) * .e17) - .e55 * .e3/((.e13) * .e17)^2)
  .e57 <- (S * (0.5 * .e54 - .e56 * .e33 * .e30) * epsilon2 - 0.5 * (.e30 * .e25/(.e13)))
  .e58 <- (.e42 * (1 + .e9 + .e10 + .e11) * sqrt(.e13))
  .e59 <- (0.5 * ((1 - .e4/(.e13)) * .e3/.e17) + .e17)
  .e60 <- (.e59 * .e33 * .e30 * .e3/((.e13) * .e17)^2 + 0.5 * .e54)
  .e61 <- (S * .e60 * epsilon2 - 0.5 * (.e30 * .e25/(.e13)))
  .e62 <- (.e34 * .e31 * .e5/.e18 + S * .e31 * .e26 * epsilon3)
  .e63 <- (.e42 * (1 + .e9 + .e10 + .e11) * (.e14) * sqrt(.e14))
  .e64 <- (S * .e31 * .e26 * epsilon3/(.e14)^2)
  .e65 <- (0.5 * ((1 - .e5/(.e14)) * .e6/.e18) + .e18)
  .e66 <- (1/((.e14) * .e18) - .e65 * .e5/((.e14) * .e18)^2)
  .e67 <- (S * (0.5 * .e64 - .e66 * .e34 * .e31) * epsilon3 - 0.5 * (.e31 * .e26/(.e14)))
  .e68 <- (.e42 * (1 + .e9 + .e10 + .e11) * sqrt(.e14))
  .e69 <- (0.5 * ((1 - .e6/(.e14)) * .e5/.e18) + .e18)
  .e70 <- (.e69 * .e34 * .e31 * .e5/((.e14) * .e18)^2 + 0.5 * .e64)
  .e71 <- (S * .e70 * epsilon3 - 0.5 * (.e31 * .e26/(.e14)))
  .e72 <- (.e35 * .e32 * .e7/.e19 + S * .e32 * .e27 * epsilon4)
  .e73 <- (.e42 * (.e15) * sqrt(.e15))
  .e74 <- (S * .e32 * .e27 * epsilon4/(.e15)^2)
  .e75 <- (0.5 * ((1 - .e7/(.e15)) * .e8/.e19) + .e19)
  .e76 <- (1/((.e15) * .e19) - .e75 * .e7/((.e15) * .e19)^2)
  .e77 <- (S * (0.5 * .e74 - .e76 * .e35 * .e32) * epsilon4 - 0.5 * (.e32 * .e27/(.e15)))
  .e78 <- (0.5 * ((1 - .e8/(.e15)) * .e7/.e19) + .e19)
  .e79 <- (.e78 * .e35 * .e32 * .e7/((.e15) * .e19)^2 + 0.5 * .e74)
  .e80 <- (S * .e79 * epsilon4 - 0.5 * (.e32 * .e27/(.e15)))
  .e81 <- (2 * (.e29 * .e24/sqrt(.e12)) - .e42)
  .e82 <- (.e42 * (1 + .e9 + .e10 + .e11))
  .e83 <- (2 * (.e30 * .e25/sqrt(.e13)) - .e42)
  .e84 <- (2 * (.e31 * .e26/sqrt(.e14)) - .e42)
  .e85 <- ((2 - 2 * (.e9/(1 + .e9 + .e10 + .e11)))/.e82 - 2 * (.e81 * .e9/.e82^2))
  .e86 <- (2 * (.e83/.e82^2) + 2/(.e42 * (1 + .e9 + .e10 + .e11)^2))
  .e87 <- (2 * (.e84/.e82^2) + 2/(.e42 * (1 + .e9 + .e10 + .e11)^2))
  .e88 <- ((2 - 2 * (.e10/(1 + .e9 + .e10 + .e11)))/.e82 - 2 * (.e83 * .e10/.e82^2))
  .e89 <- ((2 - 2 * (.e11/(1 + .e9 + .e10 + .e11)))/.e82 - 2 * (.e84 * .e11/.e82^2))
  .e90 <- (2 * ((1 + .e9 + .e10 + .e11) * .e81/.e82^2) + 2/.e82)
  .e91 <- (2 * ((1 + .e9 + .e10 + .e11) * .e83/.e82^2) + 2/.e82)
  .e92 <- (2 * ((1 + .e9 + .e10 + .e11) * .e84/.e82^2) + 2/.e82)
  .e93 <- (S * (.e28 * .e1/.e16 + S * .e24 * epsilon1) * epsilon1/(.e12) - .e24)
  hessll <- matrix(nrow = (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar),
    ncol = (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar))
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * wHvar * 2 * (((.e29 * .e93 + S *
    .e28 * (.e29 * .e1/.e2 + .e29) * .e1 * epsilon1/((.e12) * .e16))/.e43 - 2 *
    (.e39^2 * .e9/.e43^2)) * .e9), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(Xvar * wHvar * 2 *
    (S * ((.e45 * .e28 * .e29 + (S * ((0.5 * .e93 - 0.5 * .e24) * .e29/(.e12) -
      S * .e45 * .e28 * (.e29 * .e1/.e2 + .e29) * epsilon1) * epsilon1 - 0.5 *
      (.e39/(.e12)))/(.e12))/.e48 - 2 * (.e39 * .e9 * .e47/(.e48^2 * (.e12)))) *
      .e1 * .e9), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * .e93 - 0.5 * .e24) * .e29/(.e12) + S * .e49 *
    .e28 * (.e29 * .e1/.e2 + .e29) * .e1 * epsilon1/((.e12) * .e16)^2) * epsilon1 -
    0.5 * (.e39/(.e12)))/(.e12) - .e49 * .e28 * .e29 * .e1/((.e12) * .e16)^2)/.e48 -
    2 * (.e39 * .e9 * .e51/(.e48^2 * (.e12)))) * .e2 * .e9), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e39 * .e52 * (.e13) * .e9 * .e10 * sqrt(.e13)/(.e53^2 *
    (.e12) * sqrt(.e12))))), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e39 * .e3 * .e9 * .e10 *
    .e57 * sqrt(.e13)/(.e58^2 * (.e12) * sqrt(.e12))))), uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e39 * .e4 * .e9 *
    .e10 * .e61 * sqrt(.e13)/(.e58^2 * (.e12) * sqrt(.e12))))), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e39 * .e62 *
    (.e14) * .e9 * .e11 * sqrt(.e14)/(.e63^2 * (.e12) * sqrt(.e12))))), Xvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e39 * .e5 *
    .e9 * .e11 * .e67 * sqrt(.e14)/(.e68^2 * (.e12) * sqrt(.e12))))), uHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e39 * .e6 *
    .e9 * .e11 * .e71 * sqrt(.e14)/(.e68^2 * (.e12) * sqrt(.e12))))), vHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e40 * .e39 *
    .e72 * (.e15) * .e9 * sqrt(.e15)/(.e73^2 * (1 + .e9 + .e10 + .e11) * (.e12) *
    sqrt(.e12))))), Xvar)
  hessll[1:nXvar, (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e40 * .e39 *
    .e7 * .e9 * .e77 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 +
    .e11) * (.e12) * sqrt(.e12))))), uHvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e40 * .e39 *
    .e8 * .e9 * .e80 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 +
    .e11) * (.e12) * sqrt(.e12))))), vHvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * S * .e85 * .e39 *
    .e9/((.e12) * sqrt(.e12)), Zvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e86 * .e39 * .e9 * .e10/((.e12) * sqrt(.e12)))), Zvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar *
    (-(S * .e87 * .e39 * .e9 * .e11/((.e12) * sqrt(.e12)))), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e1 * (S * (0.5 * (S * .e29 * (S * (0.5 * (S * .e24 * epsilon1/(.e12)^2) -
    .e45 * .e28) * epsilon1 - 2 * (.e24/(.e12))) * epsilon1/(.e12)^2) - (0.5 *
    (S^2 * .e45 * .e29 * epsilon1^2/(.e12)^2) - (((0.5 * (.e1/(.e12)) + 1 - 0.5 *
    (0.5 * (1 - .e1/(.e12)) + .e1/(.e12))) * (1 - .e1/(.e12)) * .e2/.e16 + (2 -
    2 * (.e44^2 * .e1 * (.e12)/((.e12) * .e16)^2)) * .e16)/((.e12) * .e16)^2 +
    S^2 * .e45^2 * .e1 * epsilon1^2/((.e12) * .e16)) * .e29) * .e28) * epsilon1 -
    0.5 * ((S * (0.5 * .e46 - .e45 * .e28 * .e29) * epsilon1 - .e29 * .e24/(.e12))/(.e12))) +
    S * (0.5 * .e46 - .e45 * .e28 * .e29) * epsilon1 - 0.5 * (.e29 * .e24/(.e12)))/.e48 -
    (0.5 * (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e12)) + 2 * (.e9 * .e47)) *
      .e1 * .e47/.e48^2) * .e1 * .e9), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 * ((1 - .e1/(.e12)) *
    .e2) - S^2 * .e49 * .e45 * .e1 * epsilon1^2)/(.e12) + 0.5 * ((.e1/(.e12) -
    1) * .e2/(.e12) + 1 - 0.5 * ((1 - .e1/(.e12)) * (1 - .e2/(.e12))))) * .e29/.e16 +
    0.5 * (S^2 * .e49 * .e29 * epsilon1^2/(.e12)^2)) * .e1 + .e49 * (1 - 2 *
    (.e44 * .e1 * (.e12) * .e16/((.e12) * .e16)^2)) * .e29) * .e28/((.e12) *
    .e16)^2 + 0.5 * (S * .e29 * (S * (0.5 * (S * .e24 * epsilon1/(.e12)^2) -
    .e45 * .e28) * epsilon1 - 2 * (.e24/(.e12))) * epsilon1/(.e12)^2)) * epsilon1 -
    0.5 * ((S * (0.5 * .e46 - .e45 * .e28 * .e29) * epsilon1 - .e29 * .e24/(.e12))/(.e12)))/.e48 -
    (0.5 * (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e12)) + 2 * (.e9 * .e47)) *
      .e51/.e48^2) * .e1 * .e2 * .e9), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e52 * .e1 *
    (.e13) * .e9 * .e10 * .e47 * sqrt(.e13)/(.e53^2 * sqrt(.e12))))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e1 *
    .e3 * .e9 * .e10 * .e47 * .e57 * sqrt(.e13)/(.e58^2 * sqrt(.e12))))), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e1 *
    .e4 * .e9 * .e10 * .e61 * .e47 * sqrt(.e13)/(.e58^2 * sqrt(.e12))))), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (S * .e62 * .e1 * (.e14) * .e9 * .e11 * .e47 * sqrt(.e14)/(.e63^2 *
      sqrt(.e12))))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e1 * .e5 * .e9 * .e11 * .e47 * .e67 * sqrt(.e14)/(.e68^2 * sqrt(.e12))))),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e1 * .e6 * .e9 * .e11 * .e71 * .e47 * sqrt(.e14)/(.e68^2 * sqrt(.e12))))),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (S * .e40 * .e72 * .e1 * (.e15) * .e9 * .e47 * sqrt(.e15)/(.e73^2 *
      (1 + .e9 + .e10 + .e11) * sqrt(.e12))))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e40 * .e1 * .e7 * .e9 * .e47 * .e77 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 *
      (1 + .e9 + .e10 + .e11) * sqrt(.e12))))), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e40 * .e1 * .e8 * .e9 * .e80 * .e47 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 *
      (1 + .e9 + .e10 + .e11) * sqrt(.e12))))), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(uHvar *
    wHvar * .e85 * .e1 * .e9 * .e47/sqrt(.e12), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar *
    wHvar * (-(.e86 * .e1 * .e9 * .e10 * .e47/sqrt(.e12))), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar *
    wHvar * (-(.e87 * .e1 * .e9 * .e11 * .e47/sqrt(.e12))), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e2/(.e12)) - 0.5 * (0.5 * (1 - .e2/(.e12)) + .e2/(.e12))) * (1 - .e2/(.e12)) +
    S^2 * .e49^2 * .e1 * .e2 * epsilon1^2/(((.e12) * .e16)^2 * (.e12))) * .e29 *
    .e1/.e16 + ((0.5 * (S^2 * .e29 * epsilon1^2/(.e12)^2) - 2 * (.e49 * .e29 *
    (.e12) * .e16/((.e12) * .e16)^2)) * .e2 + .e29) * .e49) * .e28 * .e1/((.e12) *
    .e16)^2 + S * (0.5 * (.e2 * (S * (.e49 * .e28 * .e1/((.e12) * .e16)^2 + 0.5 *
    (S * .e24 * epsilon1/(.e12)^2)) * epsilon1 - 2 * (.e24/(.e12)))) + 0.5 *
    .e24) * .e29 * epsilon1/(.e12)^2) * epsilon1 - (0.5 * (.e29 * .e24) + 0.5 *
    (.e2 * (S * .e50 * epsilon1 - .e29 * .e24/(.e12))))/(.e12))/.e48 - (0.5 *
    (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e12)) + 2 * (.e9 * .e51)) * .e2 * .e51/.e48^2) *
    .e2 * .e9), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (S * .e52 * (.e13) * .e2 * .e9 * .e10 * .e51 * sqrt(.e13)/(.e53^2 *
      sqrt(.e12))))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (.e3 * .e2 * .e9 * .e10 * .e51 * .e57 * sqrt(.e13)/(.e58^2 * sqrt(.e12))))),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e2 * .e4 * .e9 * .e10 * .e51 * .e61 * sqrt(.e13)/(.e58^2 *
    sqrt(.e12))))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (S * .e62 * (.e14) * .e2 * .e9 * .e11 * .e51 * sqrt(.e14)/(.e63^2 *
    sqrt(.e12))))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e5 * .e2 * .e9 * .e11 * .e51 * .e67 * sqrt(.e14)/(.e68^2 *
    sqrt(.e12))))), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e2 * .e6 * .e9 * .e11 * .e51 * .e71 * sqrt(.e14)/(.e68^2 *
    sqrt(.e12))))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (S * .e40 * .e72 * (.e15) * .e2 * .e9 * .e51 * sqrt(.e15)/(.e73^2 *
    (1 + .e9 + .e10 + .e11) * sqrt(.e12))))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e40 * .e7 * .e2 * .e9 * .e51 * .e77 * sqrt(.e15)/((.e42 *
    sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) * sqrt(.e12))))), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e40 * .e2 * .e8 * .e9 * .e51 * .e80 * sqrt(.e15)/((.e42 *
    sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) * sqrt(.e12))))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(vHvar *
    wHvar * .e85 * .e2 * .e9 * .e51/sqrt(.e12), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar *
    wHvar * (-(.e86 * .e2 * .e9 * .e10 * .e51/sqrt(.e12))), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 *
    nZHvar)] <- crossprod(vHvar * wHvar * (-(.e87 * .e2 * .e9 * .e11 * .e51/sqrt(.e12))),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (((.e30 * (S * (.e33 * .e3/.e17 + S * .e25 * epsilon2) * epsilon2/(.e13) -
    .e25) + S * .e33 * (.e30 * .e3/.e4 + .e30) * .e3 * epsilon2/((.e13) * .e17))/.e53 -
    2 * (.e52^2 * .e10/.e53^2)) * .e10), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * ((.e56 * .e33 * .e30 + (S * ((0.5 * (S * (.e33 * .e3/.e17 +
    S * .e25 * epsilon2) * epsilon2/(.e13) - .e25) - 0.5 * .e25) * .e30/(.e13) -
    S * .e56 * .e33 * (.e30 * .e3/.e4 + .e30) * epsilon2) * epsilon2 - 0.5 *
    (.e52/(.e13)))/(.e13))/.e58 - 2 * (.e52 * .e10 * .e57/(.e58^2 * (.e13)))) *
    .e3 * .e10), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * (S * (.e33 * .e3/.e17 + S * .e25 * epsilon2) *
    epsilon2/(.e13) - .e25) - 0.5 * .e25) * .e30/(.e13) + S * .e59 * .e33 * (.e30 *
    .e3/.e4 + .e30) * .e3 * epsilon2/((.e13) * .e17)^2) * epsilon2 - 0.5 * (.e52/(.e13)))/(.e13) -
    .e59 * .e33 * .e30 * .e3/((.e13) * .e17)^2)/.e58 - 2 * (.e52 * .e10 * .e61/(.e58^2 *
    (.e13)))) * .e4 * .e10), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e52 * .e62 * (.e14) * .e10 * .e11 * sqrt(.e14)/(.e63^2 *
    (.e13) * sqrt(.e13))))), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e52 * .e5 * .e10 * .e11 * .e67 * sqrt(.e14)/(.e68^2 *
    (.e13) * sqrt(.e13))))), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e52 * .e6 * .e10 * .e11 * .e71 * sqrt(.e14)/(.e68^2 *
    (.e13) * sqrt(.e13))))), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e40 * .e52 * .e72 * (.e15) * .e10 * sqrt(.e15)/(.e73^2 *
    (1 + .e9 + .e10 + .e11) * (.e13) * sqrt(.e13))))), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e40 * .e52 * .e7 * .e10 * .e77 * sqrt(.e15)/((.e42 *
    sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) * (.e13) * sqrt(.e13))))), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e40 * .e52 * .e8 * .e10 * .e80 * sqrt(.e15)/((.e42 *
    sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) * (.e13) * sqrt(.e13))))), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    nZHvar)] <- crossprod(Xvar * wHvar * (-(S * (2 * (.e81/.e82^2) + 2/(.e42 *
    (1 + .e9 + .e10 + .e11)^2)) * .e52 * .e9 * .e10/((.e13) * sqrt(.e13)))),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * S * .e88 * .e52 *
    .e10/((.e13) * sqrt(.e13)), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S * .e87 * .e52 *
    .e10 * .e11/((.e13) * sqrt(.e13)))), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e3 * (S * (0.5 * (S * .e30 * (S * (0.5 * (S * .e25 * epsilon2/(.e13)^2) -
    .e56 * .e33) * epsilon2 - 2 * (.e25/(.e13))) * epsilon2/(.e13)^2) - (0.5 *
    (S^2 * .e56 * .e30 * epsilon2^2/(.e13)^2) - (((0.5 * (.e3/(.e13)) + 1 - 0.5 *
    (0.5 * (1 - .e3/(.e13)) + .e3/(.e13))) * (1 - .e3/(.e13)) * .e4/.e17 + (2 -
    2 * (.e55^2 * .e3 * (.e13)/((.e13) * .e17)^2)) * .e17)/((.e13) * .e17)^2 +
    S^2 * .e56^2 * .e3 * epsilon2^2/((.e13) * .e17)) * .e30) * .e33) * epsilon2 -
    0.5 * ((S * (0.5 * .e54 - .e56 * .e33 * .e30) * epsilon2 - .e30 * .e25/(.e13))/(.e13))) +
    S * (0.5 * .e54 - .e56 * .e33 * .e30) * epsilon2 - 0.5 * (.e30 * .e25/(.e13)))/.e58 -
    (0.5 * (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e13)) + 2 * (.e10 * .e57)) *
      .e3 * .e57/.e58^2) * .e3 * .e10), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((S * (((((0.5 * ((1 - .e3/(.e13)) * .e4) - S^2 * .e59 * .e56 *
    .e3 * epsilon2^2)/(.e13) + 0.5 * ((.e3/(.e13) - 1) * .e4/(.e13) + 1 - 0.5 *
    ((1 - .e3/(.e13)) * (1 - .e4/(.e13))))) * .e30/.e17 + 0.5 * (S^2 * .e59 *
    .e30 * epsilon2^2/(.e13)^2)) * .e3 + .e59 * (1 - 2 * (.e55 * .e3 * (.e13) *
    .e17/((.e13) * .e17)^2)) * .e30) * .e33/((.e13) * .e17)^2 + 0.5 * (S * .e30 *
    (S * (0.5 * (S * .e25 * epsilon2/(.e13)^2) - .e56 * .e33) * epsilon2 - 2 *
      (.e25/(.e13))) * epsilon2/(.e13)^2)) * epsilon2 - 0.5 * ((S * (0.5 *
    .e54 - .e56 * .e33 * .e30) * epsilon2 - .e30 * .e25/(.e13))/(.e13)))/.e58 -
    (0.5 * (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e13)) + 2 * (.e10 * .e57)) *
      .e61/.e58^2) * .e3 * .e4 * .e10), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e62 * .e3 * (.e14) *
    .e10 * .e11 * .e57 * sqrt(.e14)/(.e63^2 * sqrt(.e13))))), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e3 * .e5 * .e10 * .e11 *
    .e57 * .e67 * sqrt(.e14)/(.e68^2 * sqrt(.e13))))), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e3 * .e6 * .e10 * .e11 *
    .e71 * .e57 * sqrt(.e14)/(.e68^2 * sqrt(.e13))))), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e40 * .e72 * .e3 *
    (.e15) * .e10 * .e57 * sqrt(.e15)/(.e73^2 * (1 + .e9 + .e10 + .e11) * sqrt(.e13))))),
    Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e40 * .e3 * .e7 * .e10 *
    .e57 * .e77 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) *
    sqrt(.e13))))), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e40 * .e3 * .e8 * .e10 *
    .e80 * .e57 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) *
    sqrt(.e13))))), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 *
      nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-((2 * (.e81/.e82^2) +
    2/(.e42 * (1 + .e9 + .e10 + .e11)^2)) * .e3 * .e9 * .e10 * .e57/sqrt(.e13))),
    Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar +
      4 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * .e88 * .e3 *
    .e10 * .e57/sqrt(.e13), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 *
      nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e87 *
    .e3 * .e10 * .e11 * .e57/sqrt(.e13))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 * (.e4/(.e13)) -
    0.5 * (0.5 * (1 - .e4/(.e13)) + .e4/(.e13))) * (1 - .e4/(.e13)) + S^2 * .e59^2 *
    .e3 * .e4 * epsilon2^2/(((.e13) * .e17)^2 * (.e13))) * .e30 * .e3/.e17 +
    ((0.5 * (S^2 * .e30 * epsilon2^2/(.e13)^2) - 2 * (.e59 * .e30 * (.e13) *
      .e17/((.e13) * .e17)^2)) * .e4 + .e30) * .e59) * .e33 * .e3/((.e13) *
    .e17)^2 + S * (0.5 * (.e4 * (S * (.e59 * .e33 * .e3/((.e13) * .e17)^2 + 0.5 *
    (S * .e25 * epsilon2/(.e13)^2)) * epsilon2 - 2 * (.e25/(.e13)))) + 0.5 *
    .e25) * .e30 * epsilon2/(.e13)^2) * epsilon2 - (0.5 * (.e30 * .e25) + 0.5 *
    (.e4 * (S * .e60 * epsilon2 - .e30 * .e25/(.e13))))/(.e13))/.e58 - (0.5 *
    (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e13)) + 2 * (.e10 * .e61)) * .e4 *
    .e61/.e58^2) * .e4 * .e10), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e62 * (.e14) * .e4 *
    .e10 * .e11 * .e61 * sqrt(.e14)/(.e63^2 * sqrt(.e13))))), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e5 * .e4 * .e10 * .e11 *
    .e61 * .e67 * sqrt(.e14)/(.e68^2 * sqrt(.e13))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e4 * .e6 * .e10 * .e11 *
    .e61 * .e71 * sqrt(.e14)/(.e68^2 * sqrt(.e13))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e40 * .e72 * (.e15) *
    .e4 * .e10 * .e61 * sqrt(.e15)/(.e73^2 * (1 + .e9 + .e10 + .e11) * sqrt(.e13))))),
    Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e40 * .e7 * .e4 * .e10 *
    .e61 * .e77 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) *
    sqrt(.e13))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e40 * .e4 * .e8 * .e10 *
    .e61 * .e80 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 + .e11) *
    sqrt(.e13))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-((2 * (.e81/.e82^2) +
    2/(.e42 * (1 + .e9 + .e10 + .e11)^2)) * .e4 * .e9 * .e10 * .e61/sqrt(.e13))),
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * .e88 *
    .e4 * .e10 * .e61/sqrt(.e13), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e87 *
    .e4 * .e10 * .e11 * .e61/sqrt(.e13))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (((.e31 * (S * (.e34 *
    .e5/.e18 + S * .e26 * epsilon3) * epsilon3/(.e14) - .e26) + S * .e34 * (.e31 *
    .e5/.e6 + .e31) * .e5 * epsilon3/((.e14) * .e18))/.e63 - 2 * (.e62^2 * .e11/.e63^2)) *
    .e11), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * ((.e66 * .e34 *
    .e31 + (S * ((0.5 * (S * (.e34 * .e5/.e18 + S * .e26 * epsilon3) * epsilon3/(.e14) -
    .e26) - 0.5 * .e26) * .e31/(.e14) - S * .e66 * .e34 * (.e31 * .e5/.e6 + .e31) *
    epsilon3) * epsilon3 - 0.5 * (.e62/(.e14)))/(.e14))/.e68 - 2 * (.e62 * .e11 *
    .e67/(.e68^2 * (.e14)))) * .e5 * .e11), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * (((S * ((0.5 *
    (S * (.e34 * .e5/.e18 + S * .e26 * epsilon3) * epsilon3/(.e14) - .e26) -
    0.5 * .e26) * .e31/(.e14) + S * .e69 * .e34 * (.e31 * .e5/.e6 + .e31) * .e5 *
    epsilon3/((.e14) * .e18)^2) * epsilon3 - 0.5 * (.e62/(.e14)))/(.e14) - .e69 *
    .e34 * .e31 * .e5/((.e14) * .e18)^2)/.e68 - 2 * (.e62 * .e11 * .e71/(.e68^2 *
    (.e14)))) * .e6 * .e11), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e40 * .e62 *
    .e72 * (.e15) * .e11 * sqrt(.e15)/(.e73^2 * (1 + .e9 + .e10 + .e11) * (.e14) *
    sqrt(.e14))))), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e40 * .e62 *
    .e7 * .e11 * .e77 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 +
    .e11) * (.e14) * sqrt(.e14))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e40 * .e62 *
    .e8 * .e11 * .e80 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 + .e10 +
    .e11) * (.e14) * sqrt(.e14))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * (-(S * (2 *
    (.e81/.e82^2) + 2/(.e42 * (1 + .e9 + .e10 + .e11)^2)) * .e62 * .e9 * .e11/((.e14) *
    sqrt(.e14)))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e86 * .e62 * .e10 * .e11/((.e14) * sqrt(.e14)))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar *
    S * .e89 * .e62 * .e11/((.e14) * sqrt(.e14)), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((.e5 * (S * (0.5 *
    (S * .e31 * (S * (0.5 * (S * .e26 * epsilon3/(.e14)^2) - .e66 * .e34) * epsilon3 -
      2 * (.e26/(.e14))) * epsilon3/(.e14)^2) - (0.5 * (S^2 * .e66 * .e31 *
    epsilon3^2/(.e14)^2) - (((0.5 * (.e5/(.e14)) + 1 - 0.5 * (0.5 * (1 - .e5/(.e14)) +
    .e5/(.e14))) * (1 - .e5/(.e14)) * .e6/.e18 + (2 - 2 * (.e65^2 * .e5 * (.e14)/((.e14) *
    .e18)^2)) * .e18)/((.e14) * .e18)^2 + S^2 * .e66^2 * .e5 * epsilon3^2/((.e14) *
    .e18)) * .e31) * .e34) * epsilon3 - 0.5 * ((S * (0.5 * .e64 - .e66 * .e34 *
    .e31) * epsilon3 - .e31 * .e26/(.e14))/(.e14))) + S * (0.5 * .e64 - .e66 *
    .e34 * .e31) * epsilon3 - 0.5 * (.e31 * .e26/(.e14)))/.e68 - (0.5 * (.e42 *
    (1 + .e9 + .e10 + .e11)/sqrt(.e14)) + 2 * (.e11 * .e67)) * .e5 * .e67/.e68^2) *
    .e5 * .e11), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 *
    ((1 - .e5/(.e14)) * .e6) - S^2 * .e69 * .e66 * .e5 * epsilon3^2)/(.e14) +
    0.5 * ((.e5/(.e14) - 1) * .e6/(.e14) + 1 - 0.5 * ((1 - .e5/(.e14)) * (1 -
      .e6/(.e14))))) * .e31/.e18 + 0.5 * (S^2 * .e69 * .e31 * epsilon3^2/(.e14)^2)) *
    .e5 + .e69 * (1 - 2 * (.e65 * .e5 * (.e14) * .e18/((.e14) * .e18)^2)) * .e31) *
    .e34/((.e14) * .e18)^2 + 0.5 * (S * .e31 * (S * (0.5 * (S * .e26 * epsilon3/(.e14)^2) -
    .e66 * .e34) * epsilon3 - 2 * (.e26/(.e14))) * epsilon3/(.e14)^2)) * epsilon3 -
    0.5 * ((S * (0.5 * .e64 - .e66 * .e34 * .e31) * epsilon3 - .e31 * .e26/(.e14))/(.e14)))/.e68 -
    (0.5 * (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e14)) + 2 * (.e11 * .e67)) *
      .e71/.e68^2) * .e5 * .e6 * .e11), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e40 * .e72 *
    .e5 * (.e15) * .e11 * .e67 * sqrt(.e15)/(.e73^2 * (1 + .e9 + .e10 + .e11) *
    sqrt(.e14))))), Xvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e40 * .e5 *
    .e7 * .e11 * .e67 * .e77 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 +
    .e10 + .e11) * sqrt(.e14))))), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e40 * .e5 *
    .e8 * .e11 * .e80 * .e67 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 +
    .e10 + .e11) * sqrt(.e14))))), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-((2 * (.e81/.e82^2) +
    2/(.e42 * (1 + .e9 + .e10 + .e11)^2)) * .e5 * .e9 * .e11 * .e67/sqrt(.e14))),
    Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e86 *
    .e5 * .e10 * .e11 * .e67/sqrt(.e14))), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar * wHvar *
    .e89 * .e5 * .e11 * .e67/sqrt(.e14), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e6/(.e14)) - 0.5 * (0.5 * (1 - .e6/(.e14)) + .e6/(.e14))) * (1 - .e6/(.e14)) +
    S^2 * .e69^2 * .e5 * .e6 * epsilon3^2/(((.e14) * .e18)^2 * (.e14))) * .e31 *
    .e5/.e18 + ((0.5 * (S^2 * .e31 * epsilon3^2/(.e14)^2) - 2 * (.e69 * .e31 *
    (.e14) * .e18/((.e14) * .e18)^2)) * .e6 + .e31) * .e69) * .e34 * .e5/((.e14) *
    .e18)^2 + S * (0.5 * (.e6 * (S * (.e69 * .e34 * .e5/((.e14) * .e18)^2 + 0.5 *
    (S * .e26 * epsilon3/(.e14)^2)) * epsilon3 - 2 * (.e26/(.e14)))) + 0.5 *
    .e26) * .e31 * epsilon3/(.e14)^2) * epsilon3 - (0.5 * (.e31 * .e26) + 0.5 *
    (.e6 * (S * .e70 * epsilon3 - .e31 * .e26/(.e14))))/(.e14))/.e68 - (0.5 *
    (.e42 * (1 + .e9 + .e10 + .e11)/sqrt(.e14)) + 2 * (.e11 * .e71)) * .e6 *
    .e71/.e68^2) * .e6 * .e11), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e40 * .e72 *
    (.e15) * .e6 * .e11 * .e71 * sqrt(.e15)/(.e73^2 * (1 + .e9 + .e10 + .e11) *
    sqrt(.e14))))), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e40 * .e7 *
    .e6 * .e11 * .e71 * .e77 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 +
    .e10 + .e11) * sqrt(.e14))))), uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e40 * .e6 *
    .e8 * .e11 * .e71 * .e80 * sqrt(.e15)/((.e42 * sqrt(.e15))^2 * (1 + .e9 +
    .e10 + .e11) * sqrt(.e14))))), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-((2 * (.e81/.e82^2) +
    2/(.e42 * (1 + .e9 + .e10 + .e11)^2)) * .e6 * .e9 * .e11 * .e71/sqrt(.e14))),
    Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e86 *
    .e6 * .e10 * .e11 * .e71/sqrt(.e14))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(vHvar * wHvar *
    .e89 * .e6 * .e11 * .e71/sqrt(.e14), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (((.e32 * (S * (.e35 *
    .e7/.e19 + S * .e27 * epsilon4) * epsilon4/(.e15) - .e27) + S * .e35 * (.e32 *
    .e7/.e8 + .e32) * .e7 * epsilon4/((.e15) * .e19))/.e73 - 2 * (.e40 * .e72^2/.e73^2)) *
    .e40), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * ((.e76 * .e35 *
    .e32 + (S * ((0.5 * (S * (.e35 * .e7/.e19 + S * .e27 * epsilon4) * epsilon4/(.e15) -
    .e27) - 0.5 * .e27) * .e32/(.e15) - S * .e76 * .e35 * (.e32 * .e7/.e8 + .e32) *
    epsilon4) * epsilon4 - 0.5 * (.e72/(.e15)))/(.e15))/(.e42 * sqrt(.e15)) -
    2 * (.e40 * .e72 * .e77/((.e42 * sqrt(.e15))^2 * (.e15)))) * .e40 * .e7),
    uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * (((S * ((0.5 *
    (S * (.e35 * .e7/.e19 + S * .e27 * epsilon4) * epsilon4/(.e15) - .e27) -
    0.5 * .e27) * .e32/(.e15) + S * .e78 * .e35 * (.e32 * .e7/.e8 + .e32) * .e7 *
    epsilon4/((.e15) * .e19)^2) * epsilon4 - 0.5 * (.e72/(.e15)))/(.e15) - .e78 *
    .e35 * .e32 * .e7/((.e15) * .e19)^2)/(.e42 * sqrt(.e15)) - 2 * (.e40 * .e72 *
    .e80/((.e42 * sqrt(.e15))^2 * (.e15)))) * .e40 * .e8), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * (-(S * .e40 *
    .e90 * .e72 * .e9/((.e15) * sqrt(.e15)))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e40 * .e91 * .e72 * .e10/((.e15) * sqrt(.e15)))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar *
    (-(S * .e40 * .e92 * .e72 * .e11/((.e15) * sqrt(.e15)))), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((.e7 * (S * (0.5 *
    (S * .e32 * (S * (0.5 * (S * .e27 * epsilon4/(.e15)^2) - .e76 * .e35) * epsilon4 -
      2 * (.e27/(.e15))) * epsilon4/(.e15)^2) - (0.5 * (S^2 * .e76 * .e32 *
    epsilon4^2/(.e15)^2) - (((0.5 * (.e7/(.e15)) + 1 - 0.5 * (0.5 * (1 - .e7/(.e15)) +
    .e7/(.e15))) * (1 - .e7/(.e15)) * .e8/.e19 + (2 - 2 * (.e75^2 * .e7 * (.e15)/((.e15) *
    .e19)^2)) * .e19)/((.e15) * .e19)^2 + S^2 * .e76^2 * .e7 * epsilon4^2/((.e15) *
    .e19)) * .e32) * .e35) * epsilon4 - 0.5 * ((S * (0.5 * .e74 - .e76 * .e35 *
    .e32) * epsilon4 - .e32 * .e27/(.e15))/(.e15))) + S * (0.5 * .e74 - .e76 *
    .e35 * .e32) * epsilon4 - 0.5 * (.e32 * .e27/(.e15)))/(.e42 * sqrt(.e15)) -
    (0.5 * (.e42/sqrt(.e15)) + 2 * (.e40 * .e77)) * .e7 * .e77/(.e42 * sqrt(.e15))^2) *
    .e40 * .e7), uHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 *
    ((1 - .e7/(.e15)) * .e8) - S^2 * .e78 * .e76 * .e7 * epsilon4^2)/(.e15) +
    0.5 * ((.e7/(.e15) - 1) * .e8/(.e15) + 1 - 0.5 * ((1 - .e7/(.e15)) * (1 -
      .e8/(.e15))))) * .e32/.e19 + 0.5 * (S^2 * .e78 * .e32 * epsilon4^2/(.e15)^2)) *
    .e7 + .e78 * (1 - 2 * (.e75 * .e7 * (.e15) * .e19/((.e15) * .e19)^2)) * .e32) *
    .e35/((.e15) * .e19)^2 + 0.5 * (S * .e32 * (S * (0.5 * (S * .e27 * epsilon4/(.e15)^2) -
    .e76 * .e35) * epsilon4 - 2 * (.e27/(.e15))) * epsilon4/(.e15)^2)) * epsilon4 -
    0.5 * ((S * (0.5 * .e74 - .e76 * .e35 * .e32) * epsilon4 - .e32 * .e27/(.e15))/(.e15)))/(.e42 *
    sqrt(.e15)) - (0.5 * (.e42/sqrt(.e15)) + 2 * (.e40 * .e77)) * .e80/(.e42 *
    sqrt(.e15))^2) * .e40 * .e7 * .e8), vHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-(.e40 * .e90 *
    .e7 * .e9 * .e77/sqrt(.e15))), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e40 *
    .e91 * .e7 * .e10 * .e77/sqrt(.e15))), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar * wHvar *
    (-(.e40 * .e92 * .e7 * .e11 * .e77/sqrt(.e15))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e8/(.e15)) - 0.5 * (0.5 * (1 - .e8/(.e15)) + .e8/(.e15))) * (1 - .e8/(.e15)) +
    S^2 * .e78^2 * .e7 * .e8 * epsilon4^2/(((.e15) * .e19)^2 * (.e15))) * .e32 *
    .e7/.e19 + ((0.5 * (S^2 * .e32 * epsilon4^2/(.e15)^2) - 2 * (.e78 * .e32 *
    (.e15) * .e19/((.e15) * .e19)^2)) * .e8 + .e32) * .e78) * .e35 * .e7/((.e15) *
    .e19)^2 + S * (0.5 * (.e8 * (S * (.e78 * .e35 * .e7/((.e15) * .e19)^2 + 0.5 *
    (S * .e27 * epsilon4/(.e15)^2)) * epsilon4 - 2 * (.e27/(.e15)))) + 0.5 *
    .e27) * .e32 * epsilon4/(.e15)^2) * epsilon4 - (0.5 * (.e32 * .e27) + 0.5 *
    (.e8 * (S * .e79 * epsilon4 - .e32 * .e27/(.e15))))/(.e15))/(.e42 * sqrt(.e15)) -
    (0.5 * (.e42/sqrt(.e15)) + 2 * (.e40 * .e80)) * .e8 * .e80/(.e42 * sqrt(.e15))^2) *
    .e40 * .e8), vHvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-(.e40 * .e90 *
    .e8 * .e9 * .e80/sqrt(.e15))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e40 *
    .e91 * .e8 * .e10 * .e80/sqrt(.e15))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(vHvar * wHvar *
    (-(.e40 * .e92 * .e8 * .e11 * .e80/sqrt(.e15))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + nZHvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(Zvar * wHvar * ((1 - .e9/(1 +
    .e9 + .e10 + .e11))/.e82 - 2 * (.e29 * .e9 * .e24/(.e82^2 * sqrt(.e12)))) *
    .e81 * .e9, Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + nZHvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e81/(.e42 * (1 + .e9 + .e10 + .e11)^2) + 2 * (.e83 * .e29 *
    .e24/(.e82^2 * sqrt(.e12)))) * .e9 * .e10)), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + nZHvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e81/(.e42 * (1 + .e9 + .e10 + .e11)^2) + 2 * (.e84 * .e29 *
    .e24/(.e82^2 * sqrt(.e12)))) * .e9 * .e11)), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + 2 * nZHvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(Zvar *
    wHvar * ((1 - .e10/(1 + .e9 + .e10 + .e11))/.e82 - 2 * (.e30 * .e10 * .e25/(.e82^2 *
    sqrt(.e13)))) * .e83 * .e10, Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar + 2 * nZHvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e83/(.e42 * (1 + .e9 + .e10 + .e11)^2) + 2 * (.e84 * .e30 *
    .e25/(.e82^2 * sqrt(.e13)))) * .e10 * .e11)), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(Zvar *
    wHvar * ((1 - .e11/(1 + .e9 + .e10 + .e11))/.e82 - 2 * (.e31 * .e11 * .e26/(.e82^2 *
    sqrt(.e14)))) * .e84 * .e11, Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for lcm 4 classes halfnormal-normal distribution
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
LCM4ChnormAlgOpt <- function(start, olsParam, dataTable, S, wHvar,
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csLCMfhalfnorm4C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(cLCMhalfnormlike4C(startVal, nXvar = nXvar,
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
  cat("LCM 4 Classes Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cLCMhalfnormlike4C(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cLCMhalfnormlike4C,
    grad = cgradLCMhalfnormlike4C, hess = chessLCMhalfnormlike4C,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cLCMhalfnormlike4C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessLCMhalfnormlike4C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradLCMhalfnormlike4C(mleObj$par,
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
      mleObj$hessian <- chessLCMhalfnormlike4C(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessLCMhalfnormlike4C(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cLCMhalfnormlike4C(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradLCMhalfnormlike4C(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) initHalf = initHalf))
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for lcmcross 4 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @param level level for confidence interval
#' @noRd
cLCM4Chalfnormeff <- function(object, level) {
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
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar)]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + object$nZHvar + 1):(4 * object$nXvar +
    4 * object$nuZUvar + 4 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 2 * object$nZHvar + 1):(4 * object$nXvar +
    4 * object$nuZUvar + 4 * object$nvZVvar + 3 * object$nZHvar)]
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
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
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
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4),
    1, which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c ==
    2, Pcond_c2, ifelse(Group_c == 3, Pcond_c3, Pcond_c4)))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
  u_c4 <- mustar4 + sigmastar4 * dnorm(mustar4/sigmastar4)/pnorm(mustar4/sigmastar4)
  u_c <- ifelse(Group_c == 1, u_c1, ifelse(Group_c == 2, u_c2,
    ifelse(Group_c == 3, u_c3, u_c4)))
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
    teBC_c4 <- exp(-mustar4 + 1/2 * sigmastar4^2) * pnorm(mustar4/sigmastar4 -
      sigmastar4)/pnorm(mustar4/sigmastar4)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, ifelse(Group_c ==
      2, teBC_c2, ifelse(Group_c == 3, teBC_c3, teBC_c4)))
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    effBC_c3 <- ifelse(Group_c == 3, teBC_c3, NA)
    effBC_c4 <- ifelse(Group_c == 4, teBC_c4, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) *
      pnorm(mustar1/sigmastar1 + sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) *
      pnorm(mustar2/sigmastar2 + sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c3 <- exp(mustar3 + 1/2 * sigmastar3^2) *
      pnorm(mustar3/sigmastar3 + sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_reciprocal_c4 <- exp(mustar4 + 1/2 * sigmastar4^2) *
      pnorm(mustar4/sigmastar4 + sigmastar4)/pnorm(mustar4/sigmastar4)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      ifelse(Group_c == 2, teBC_reciprocal_c2, ifelse(Group_c ==
        3, teBC_reciprocal_c3, teBC_reciprocal_c4)))
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    ReffBC_c3 <- ifelse(Group_c == 3, teBC_reciprocal_c3,
      NA)
    ReffBC_c4 <- ifelse(Group_c == 4, teBC_reciprocal_c4,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, PosteriorProb_c3 = Pcond_c3,
      PriorProb_c3 = Probc3, u_c3 = u_c3, teBC_c3 = teBC_c3,
      teBC_reciprocal_c3 = teBC_reciprocal_c3, PosteriorProb_c4 = Pcond_c4,
      PriorProb_c4 = Probc4, u_c4 = u_c4, teBC_c4 = teBC_c4,
      teBC_reciprocal_c4 = teBC_reciprocal_c4, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, ineff_c4 = ineff_c4,
      effBC_c1 = effBC_c1, effBC_c2 = effBC_c2, effBC_c3 = effBC_c3,
      effBC_c4 = effBC_c4, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2,
      ReffBC_c3 = ReffBC_c3, ReffBC_c4 = ReffBC_c4)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
      u_c3 = u_c3, PosteriorProb_c4 = Pcond_c4, PriorProb_c4 = Probc4,
      u_c4 = u_c4, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2,
      ineff_c3 = ineff_c3, ineff_c4 = ineff_c4)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal effects for for lcmcross 4 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @noRd
cmargLCM4Chalfnorm_Eu <- function(object) {
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
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar)]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + object$nZHvar + 1):(4 * object$nXvar +
    4 * object$nuZUvar + 4 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 2 * object$nZHvar + 1):(4 * object$nXvar +
    4 * object$nuZUvar + 4 * object$nvZVvar + 3 * object$nZHvar)]
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
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
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
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4),
    1, which.max)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu1/2) * dnorm(0), ncol = 1))
  colnames(margEff_c1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu2/2) * dnorm(0), ncol = 1))
  colnames(margEff_c2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu3/2) * dnorm(0), ncol = 1))
  colnames(margEff_c3) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c3")
  margEff_c4 <- kronecker(matrix(delta4[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu4/2) * dnorm(0), ncol = 1))
  colnames(margEff_c4) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c4")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      ifelse(Group_c == 2, margEff_c2[, c], ifelse(Group_c ==
        3, margEff_c3[, c], margEff_c4[, c])))
  }
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2,
    margEff_c3, margEff_c4)
  return(margEff)
}

cmargLCM4Chalfnorm_Vu <- function(object) {
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
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar)]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + object$nZHvar + 1):(4 * object$nXvar +
    4 * object$nuZUvar + 4 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 2 * object$nZHvar + 1):(4 * object$nXvar +
    4 * object$nuZUvar + 4 * object$nvZVvar + 3 * object$nZHvar)]
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
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
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
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4),
    1, which.max)
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
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu3) * (1 - (dnorm(0)/pnorm(0))^2),
    ncol = 1))
  colnames(margEff_c3) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c3")
  margEff_c4 <- kronecker(matrix(delta4[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu4) * (1 - (dnorm(0)/pnorm(0))^2),
    ncol = 1))
  colnames(margEff_c4) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c4")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      ifelse(Group_c == 2, margEff_c2[, c], ifelse(Group_c ==
        3, margEff_c3[, c], margEff_c4[, c])))
  }
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2,
    margEff_c3, margEff_c4)
  return(margEff)
}
