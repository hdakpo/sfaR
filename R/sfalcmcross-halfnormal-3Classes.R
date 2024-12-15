################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Latent Class Stochastic Frontier Analysis                             #
# Number of Classes: 3L                                                        #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for lcm 3 classes halfnormal-normal distribution
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
cLCMhalfnormlike3C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  theta1 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar)]
  theta2 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
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
  mustar3 <- -exp(Wu3) * S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2))
  Probc3 <- 1 - Probc1 - Probc2
  ll <- log(Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3)
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for lcm 3 classes halfnormal-normal distribution
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
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
csLCMfhalfnorm3C <- function(olsObj, epsiRes, nXvar, nuZUvar,
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
    rep(0, 2 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)),
    paste0("Zv_", colnames(vHvar)), paste0("Cl1_", colnames(Zvar)),
    paste0("Cl2_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for lcm 3 classes halfnormal-normal distribution
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
cgradLCMhalfnormlike3C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  theta1 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar)]
  theta2 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
   .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wu3)
  .e6 <- exp(Wv3)
  .e7 <- .e1 + .e2
  .e8 <- .e3 + .e4
  .e9 <- .e5 + .e6
  .e10 <- sqrt(.e7)
  .e11 <- sqrt(.e8)
  .e12 <- sqrt(.e9)
  .e13 <- sqrt(.e1 * .e2/(.e7))
  .e14 <- sqrt(.e5 * .e6/(.e9))
  .e15 <- sqrt(.e3 * .e4/(.e8))
  .e16 <- dnorm(S * epsilon1/.e10)
  .e17 <- dnorm(S * epsilon2/.e11)
  .e18 <- dnorm(S * epsilon3/.e12)
  .e19 <- (S * .e1 * epsilon1/((.e7) * .e13))
  .e20 <- (S * .e3 * epsilon2/((.e8) * .e15))
  .e21 <- (S * .e5 * epsilon3/((.e9) * .e14))
  .e22 <- dnorm(-.e19)
  .e23 <- pnorm(-.e19)
  .e24 <- dnorm(-.e20)
  .e25 <- pnorm(-.e20)
  .e26 <- dnorm(-.e21)
  .e27 <- pnorm(-.e21)
  .e28 <- exp(Wz1)
  .e29 <- exp(Wz2)
  .e30 <- (1 - (.e28 + .e29)/(1 + .e28 + .e29))
  .e31 <- (2 * (.e16 * .e28 * .e23/.e10) + 2 * (.e17 * .e29 * .e25/.e11))
  .e32 <- (.e22 * .e16 * .e1/.e13 + S * .e16 * .e23 * epsilon1)
  .e33 <- (.e31/(1 + .e28 + .e29) + 2 * (.e30 * .e18 * .e27/.e12))
  .e34 <- (.e33 * (1 + .e28 + .e29) * (.e7) * .e10)
  .e35 <- (S * .e16 * .e23 * epsilon1/(.e7)^2)
  .e36 <- (0.5 * ((1 - .e1/(.e7)) * .e2/.e13) + .e13)
  .e37 <- (1/((.e7) * .e13) - .e36 * .e1/((.e7) * .e13)^2)
  .e38 <- (S * (0.5 * .e35 - .e37 * .e22 * .e16) * epsilon1 - 0.5 * (.e16 * .e23/(.e7)))
  .e39 <- (.e33 * (1 + .e28 + .e29) * .e10)
  .e40 <- (0.5 * ((1 - .e2/(.e7)) * .e1/.e13) + .e13)
  .e41 <- (.e40 * .e22 * .e16 * .e1/((.e7) * .e13)^2 + 0.5 * .e35)
  .e42 <- (S * .e41 * epsilon1 - 0.5 * (.e16 * .e23/(.e7)))
  .e43 <- (.e24 * .e17 * .e3/.e15 + S * .e17 * .e25 * epsilon2)
  .e44 <- (.e33 * (1 + .e28 + .e29) * (.e8) * .e11)
  .e45 <- (S * .e17 * .e25 * epsilon2/(.e8)^2)
  .e46 <- (0.5 * ((1 - .e3/(.e8)) * .e4/.e15) + .e15)
  .e47 <- (1/((.e8) * .e15) - .e46 * .e3/((.e8) * .e15)^2)
  .e48 <- (S * (0.5 * .e45 - .e47 * .e24 * .e17) * epsilon2 - 0.5 * (.e17 * .e25/(.e8)))
  .e49 <- (.e33 * (1 + .e28 + .e29) * .e11)
  .e50 <- (0.5 * ((1 - .e4/(.e8)) * .e3/.e15) + .e15)
  .e51 <- (.e50 * .e24 * .e17 * .e3/((.e8) * .e15)^2 + 0.5 * .e45)
  .e52 <- (S * .e51 * epsilon2 - 0.5 * (.e17 * .e25/(.e8)))
  .e53 <- (.e26 * .e18 * .e5/.e14 + S * .e18 * .e27 * epsilon3)
  .e54 <- (S * .e18 * .e27 * epsilon3/(.e9)^2)
  .e55 <- (0.5 * ((1 - .e5/(.e9)) * .e6/.e14) + .e14)
  .e56 <- (1/((.e9) * .e14) - .e55 * .e5/((.e9) * .e14)^2)
  .e57 <- (S * (0.5 * .e54 - .e56 * .e26 * .e18) * epsilon3 - 0.5 * (.e18 * .e27/(.e9)))
  .e58 <- (0.5 * ((1 - .e6/(.e9)) * .e5/.e14) + .e14)
  .e59 <- (.e58 * .e26 * .e18 * .e5/((.e9) * .e14)^2 + 0.5 * .e54)
  .e60 <- (S * .e59 * epsilon3 - 0.5 * (.e18 * .e27/(.e9)))
  gradll <- cbind(2 * (S * Xvar * .e32 * .e28/.e34), 2 * (uHvar * .e1 * .e28 *
    .e38/.e39), 2 * (vHvar * .e2 * .e28 * .e42/.e39), 2 * (S * Xvar * .e43 *
    .e29/.e44), 2 * (uHvar * .e3 * .e29 * .e48/.e49), 2 * (vHvar * .e4 * .e29 *
    .e52/.e49), 2 * (S * Xvar * .e30 * .e53/(.e33 * (.e9) * .e12)), 2 * (uHvar *
    .e30 * .e5 * .e57/(.e33 * .e12)), 2 * (vHvar * .e30 * .e6 * .e60/(.e33 *
    .e12)), Zvar * (2 * (.e16 * .e23/.e10) - .e33) * .e28/(.e33 * (1 + .e28 +
    .e29)), Zvar * (2 * (.e17 * .e25/.e11) - .e33) * .e29/(.e33 * (1 + .e28 +
    .e29)))
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for lcm 3 classes halfnormal-normal distribution
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
chessLCMhalfnormlike3C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  theta1 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar)]
  theta2 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wu3)
  .e6 <- exp(Wv3)
  .e7 <- .e1 + .e2
  .e8 <- .e3 + .e4
  .e9 <- .e5 + .e6
  .e10 <- sqrt(.e7)
  .e11 <- sqrt(.e8)
  .e12 <- sqrt(.e9)
  .e13 <- sqrt(.e1 * .e2/(.e7))
  .e14 <- sqrt(.e5 * .e6/(.e9))
  .e15 <- sqrt(.e3 * .e4/(.e8))
  .e16 <- dnorm(S * epsilon1/.e10)
  .e17 <- dnorm(S * epsilon2/.e11)
  .e18 <- dnorm(S * epsilon3/.e12)
  .e19 <- (S * .e1 * epsilon1/((.e7) * .e13))
  .e20 <- (S * .e3 * epsilon2/((.e8) * .e15))
  .e21 <- (S * .e5 * epsilon3/((.e9) * .e14))
  .e22 <- dnorm(-.e19)
  .e23 <- pnorm(-.e19)
  .e24 <- dnorm(-.e20)
  .e25 <- pnorm(-.e20)
  .e26 <- dnorm(-.e21)
  .e27 <- pnorm(-.e21)
  .e28 <- exp(Wz1)
  .e29 <- exp(Wz2)
  .e30 <- (1 - (.e28 + .e29)/(1 + .e28 + .e29))
  .e31 <- (2 * (.e16 * .e28 * .e23/.e10) + 2 * (.e17 * .e29 * .e25/.e11))
  .e32 <- (.e22 * .e16 * .e1/.e13 + S * .e16 * .e23 * epsilon1)
  .e33 <- (.e31/(1 + .e28 + .e29) + 2 * (.e30 * .e18 * .e27/.e12))
  .e34 <- (.e33 * (1 + .e28 + .e29) * (.e7) * .e10)
  .e35 <- (S * .e16 * .e23 * epsilon1/(.e7)^2)
  .e36 <- (0.5 * ((1 - .e1/(.e7)) * .e2/.e13) + .e13)
  .e37 <- (1/((.e7) * .e13) - .e36 * .e1/((.e7) * .e13)^2)
  .e38 <- (S * (0.5 * .e35 - .e37 * .e22 * .e16) * epsilon1 - 0.5 * (.e16 * .e23/(.e7)))
  .e39 <- (.e33 * (1 + .e28 + .e29) * .e10)
  .e40 <- (0.5 * ((1 - .e2/(.e7)) * .e1/.e13) + .e13)
  .e41 <- (.e40 * .e22 * .e16 * .e1/((.e7) * .e13)^2 + 0.5 * .e35)
  .e42 <- (S * .e41 * epsilon1 - 0.5 * (.e16 * .e23/(.e7)))
  .e43 <- (.e24 * .e17 * .e3/.e15 + S * .e17 * .e25 * epsilon2)
  .e44 <- (.e33 * (1 + .e28 + .e29) * (.e8) * .e11)
  .e45 <- (S * .e17 * .e25 * epsilon2/(.e8)^2)
  .e46 <- (0.5 * ((1 - .e3/(.e8)) * .e4/.e15) + .e15)
  .e47 <- (1/((.e8) * .e15) - .e46 * .e3/((.e8) * .e15)^2)
  .e48 <- (S * (0.5 * .e45 - .e47 * .e24 * .e17) * epsilon2 - 0.5 * (.e17 * .e25/(.e8)))
  .e49 <- (.e33 * (1 + .e28 + .e29) * .e11)
  .e50 <- (0.5 * ((1 - .e4/(.e8)) * .e3/.e15) + .e15)
  .e51 <- (.e50 * .e24 * .e17 * .e3/((.e8) * .e15)^2 + 0.5 * .e45)
  .e52 <- (S * .e51 * epsilon2 - 0.5 * (.e17 * .e25/(.e8)))
  .e53 <- (.e26 * .e18 * .e5/.e14 + S * .e18 * .e27 * epsilon3)
  .e54 <- (S * .e18 * .e27 * epsilon3/(.e9)^2)
  .e55 <- (0.5 * ((1 - .e5/(.e9)) * .e6/.e14) + .e14)
  .e56 <- (1/((.e9) * .e14) - .e55 * .e5/((.e9) * .e14)^2)
  .e57 <- (S * (0.5 * .e54 - .e56 * .e26 * .e18) * epsilon3 - 0.5 * (.e18 * .e27/(.e9)))
  .e58 <- (0.5 * ((1 - .e6/(.e9)) * .e5/.e14) + .e14)
  .e59 <- (.e58 * .e26 * .e18 * .e5/((.e9) * .e14)^2 + 0.5 * .e54)
  .e60 <- (S * .e59 * epsilon3 - 0.5 * (.e18 * .e27/(.e9)))
  .e61 <- (S * (.e22 * .e1/.e13 + S * .e23 * epsilon1) * epsilon1/(.e7) - .e23)
  .e62 <- ((.e33 * .e12)^2 * (1 + .e28 + .e29) * (.e7) * .e10)
  .e63 <- ((2 * (.e16 * .e23/.e10) - .e33) * .e28/(.e33 * (1 + .e28 + .e29))^2)
  .e64 <- ((2 - 2 * (.e28/(1 + .e28 + .e29)))/(.e33 * (1 + .e28 + .e29)) - 2 *
    .e63)
  .e65 <- ((2 * (.e17 * .e25/.e11) - .e33)/(.e33 * (1 + .e28 + .e29))^2)
  .e66 <- (2 * .e65 + 2/(.e33 * (1 + .e28 + .e29)^2))
  .e67 <- (S * (0.5 * (S * .e23 * epsilon1/(.e7)^2) - .e37 * .e22) * epsilon1 -
    2 * (.e23/(.e7)))
  .e68 <- ((2 * (.e17 * .e25/.e11) - .e33) * .e29/(.e33 * (1 + .e28 + .e29))^2)
  .e69 <- ((2 - 2 * (.e29/(1 + .e28 + .e29)))/(.e33 * (1 + .e28 + .e29)) - 2 *
    .e68)
  .e70 <- ((1 + .e28 + .e29) * (2 * (.e16 * .e23/.e10) - .e33)/(.e33 * (1 + .e28 +
    .e29))^2)
  .e71 <- (2 * .e70 + 2/(.e33 * (1 + .e28 + .e29)))
  .e72 <- ((1 + .e28 + .e29) * (2 * (.e17 * .e25/.e11) - .e33)/(.e33 * (1 + .e28 +
    .e29))^2)
  .e73 <- (2 * .e72 + 2/(.e33 * (1 + .e28 + .e29)))
  hessll <- matrix(nrow = (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar),
    ncol = (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar))
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * wHvar * 2 * (S^2 * ((.e16 * .e61 +
    S * .e22 * (.e16 * .e1/.e2 + .e16) * .e1 * epsilon1/((.e7) * .e13))/.e34 -
    2 * (.e32^2 * .e28/.e34^2)) * .e28), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(Xvar * wHvar * 2 *
    (S * ((.e37 * .e22 * .e16 + (S * ((0.5 * .e61 - 0.5 * .e23) * .e16/(.e7) -
      S * .e37 * .e22 * (.e16 * .e1/.e2 + .e16) * epsilon1) * epsilon1 - 0.5 *
      (.e32/(.e7)))/(.e7))/.e39 - 2 * (.e32 * .e28 * .e38/(.e39^2 * (.e7)))) *
      .e1 * .e28), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * .e61 - 0.5 * .e23) * .e16/(.e7) + S * .e40 *
    .e22 * (.e16 * .e1/.e2 + .e16) * .e1 * epsilon1/((.e7) * .e13)^2) * epsilon1 -
    0.5 * (.e32/(.e7)))/(.e7) - .e40 * .e22 * .e16 * .e1/((.e7) * .e13)^2)/.e39 -
    2 * (.e32 * .e28 * .e42/(.e39^2 * (.e7)))) * .e2 * .e28), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e32 * .e43 * (.e8) * .e28 * .e29 * .e11/(.e44^2 * (.e7) *
    .e10)))), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e32 * .e3 * .e28 * .e29 *
    .e48 * .e11/(.e49^2 * (.e7) * .e10)))), uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e32 * .e4 * .e28 *
    .e29 * .e52 * .e11/(.e49^2 * (.e7) * .e10)))), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e30 * .e32 *
    .e53 * (.e9) * .e28 * .e12/((.e33 * (.e9) * .e12)^2 * (1 + .e28 + .e29) *
    (.e7) * .e10)))), Xvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e30 * .e32 *
    .e5 * .e28 * .e57 * .e12/.e62))), uHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e30 * .e32 *
    .e6 * .e28 * .e60 * .e12/.e62))), vHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * S * .e64 * .e32 *
    .e28/((.e7) * .e10), Zvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e66 * .e32 * .e28 * .e29/((.e7) * .e10))), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e1 * (S * (0.5 * (S * .e16 * .e67 * epsilon1/(.e7)^2) - (0.5 *
    (S^2 * .e37 * .e16 * epsilon1^2/(.e7)^2) - (((0.5 * (.e1/(.e7)) + 1 - 0.5 *
    (0.5 * (1 - .e1/(.e7)) + .e1/(.e7))) * (1 - .e1/(.e7)) * .e2/.e13 + (2 -
    2 * (.e36^2 * .e1 * (.e7)/((.e7) * .e13)^2)) * .e13)/((.e7) * .e13)^2 + S^2 *
    .e37^2 * .e1 * epsilon1^2/((.e7) * .e13)) * .e16) * .e22) * epsilon1 - 0.5 *
    ((S * (0.5 * .e35 - .e37 * .e22 * .e16) * epsilon1 - .e16 * .e23/(.e7))/(.e7))) +
    S * (0.5 * .e35 - .e37 * .e22 * .e16) * epsilon1 - 0.5 * (.e16 * .e23/(.e7)))/.e39 -
    (0.5 * (.e33 * (1 + .e28 + .e29)/.e10) + 2 * (.e28 * .e38)) * .e1 * .e38/.e39^2) *
    .e1 * .e28), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 * ((1 - .e1/(.e7)) *
    .e2) - S^2 * .e40 * .e37 * .e1 * epsilon1^2)/(.e7) + 0.5 * ((.e1/(.e7) -
    1) * .e2/(.e7) + 1 - 0.5 * ((1 - .e1/(.e7)) * (1 - .e2/(.e7))))) * .e16/.e13 +
    0.5 * (S^2 * .e40 * .e16 * epsilon1^2/(.e7)^2)) * .e1 + .e40 * (1 - 2 * (.e36 *
    .e1 * (.e7) * .e13/((.e7) * .e13)^2)) * .e16) * .e22/((.e7) * .e13)^2 + 0.5 *
    (S * .e16 * .e67 * epsilon1/(.e7)^2)) * epsilon1 - 0.5 * ((S * (0.5 * .e35 -
    .e37 * .e22 * .e16) * epsilon1 - .e16 * .e23/(.e7))/(.e7)))/.e39 - (0.5 *
    (.e33 * (1 + .e28 + .e29)/.e10) + 2 * (.e28 * .e38)) * .e42/.e39^2) * .e1 *
    .e2 * .e28), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e43 * .e1 *
    (.e8) * .e28 * .e29 * .e38 * .e11/(.e44^2 * .e10)))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e1 *
    .e3 * .e28 * .e29 * .e38 * .e48 * .e11/(.e49^2 * .e10)))), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e1 *
    .e4 * .e28 * .e29 * .e52 * .e38 * .e11/(.e49^2 * .e10)))), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (S * .e30 * .e53 * .e1 * (.e9) * .e28 * .e38 * .e12/((.e33 * (.e9) *
      .e12)^2 * (1 + .e28 + .e29) * .e10)))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e30 * .e1 * .e5 * .e28 * .e38 * .e57 * .e12/((.e33 * .e12)^2 * (1 +
      .e28 + .e29) * .e10)))), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e30 * .e1 * .e6 * .e28 * .e60 * .e38 * .e12/((.e33 * .e12)^2 * (1 +
      .e28 + .e29) * .e10)))), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(uHvar *
    wHvar * .e64 * .e1 * .e28 * .e38/.e10, Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    nZHvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar *
    wHvar * (-(.e66 * .e1 * .e28 * .e29 * .e38/.e10)), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e2/(.e7)) - 0.5 * (0.5 * (1 - .e2/(.e7)) + .e2/(.e7))) * (1 - .e2/(.e7)) +
    S^2 * .e40^2 * .e1 * .e2 * epsilon1^2/(((.e7) * .e13)^2 * (.e7))) * .e16 *
    .e1/.e13 + ((0.5 * (S^2 * .e16 * epsilon1^2/(.e7)^2) - 2 * (.e40 * .e16 *
    (.e7) * .e13/((.e7) * .e13)^2)) * .e2 + .e16) * .e40) * .e22 * .e1/((.e7) *
    .e13)^2 + S * (0.5 * (.e2 * (S * (.e40 * .e22 * .e1/((.e7) * .e13)^2 + 0.5 *
    (S * .e23 * epsilon1/(.e7)^2)) * epsilon1 - 2 * (.e23/(.e7)))) + 0.5 * .e23) *
    .e16 * epsilon1/(.e7)^2) * epsilon1 - (0.5 * (.e16 * .e23) + 0.5 * (.e2 *
    (S * .e41 * epsilon1 - .e16 * .e23/(.e7))))/(.e7))/.e39 - (0.5 * (.e33 *
    (1 + .e28 + .e29)/.e10) + 2 * (.e28 * .e42)) * .e2 * .e42/.e39^2) * .e2 *
    .e28), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (S * .e43 * (.e8) * .e2 * .e28 * .e29 * .e42 * .e11/(.e44^2 * .e10)))),
    Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (.e3 * .e2 * .e28 * .e29 * .e42 * .e48 * .e11/(.e49^2 * .e10)))),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e2 * .e4 * .e28 * .e29 * .e42 * .e52 * .e11/(.e49^2 * .e10)))),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (S * .e30 * .e53 * (.e9) * .e2 * .e28 * .e42 * .e12/((.e33 *
    (.e9) * .e12)^2 * (1 + .e28 + .e29) * .e10)))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e30 * .e5 * .e2 * .e28 * .e42 * .e57 * .e12/((.e33 * .e12)^2 *
    (1 + .e28 + .e29) * .e10)))), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e30 * .e2 * .e6 * .e28 * .e42 * .e60 * .e12/((.e33 * .e12)^2 *
    (1 + .e28 + .e29) * .e10)))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(vHvar *
    wHvar * .e64 * .e2 * .e28 * .e42/.e10, Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + nZHvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar *
    wHvar * (-(.e66 * .e2 * .e28 * .e29 * .e42/.e10)), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (((.e17 * (S * (.e24 * .e3/.e15 + S * .e25 * epsilon2) * epsilon2/(.e8) -
    .e25) + S * .e24 * (.e17 * .e3/.e4 + .e17) * .e3 * epsilon2/((.e8) * .e15))/.e44 -
    2 * (.e43^2 * .e29/.e44^2)) * .e29), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * ((.e47 * .e24 * .e17 + (S * ((0.5 * (S * (.e24 * .e3/.e15 +
    S * .e25 * epsilon2) * epsilon2/(.e8) - .e25) - 0.5 * .e25) * .e17/(.e8) -
    S * .e47 * .e24 * (.e17 * .e3/.e4 + .e17) * epsilon2) * epsilon2 - 0.5 *
    (.e43/(.e8)))/(.e8))/.e49 - 2 * (.e43 * .e29 * .e48/(.e49^2 * (.e8)))) *
    .e3 * .e29), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * (S * (.e24 * .e3/.e15 + S * .e25 * epsilon2) *
    epsilon2/(.e8) - .e25) - 0.5 * .e25) * .e17/(.e8) + S * .e50 * .e24 * (.e17 *
    .e3/.e4 + .e17) * .e3 * epsilon2/((.e8) * .e15)^2) * epsilon2 - 0.5 * (.e43/(.e8)))/(.e8) -
    .e50 * .e24 * .e17 * .e3/((.e8) * .e15)^2)/.e49 - 2 * (.e43 * .e29 * .e52/(.e49^2 *
    (.e8)))) * .e4 * .e29), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e30 * .e43 * .e53 * (.e9) * .e29 * .e12/((.e33 * (.e9) *
    .e12)^2 * (1 + .e28 + .e29) * (.e8) * .e11)))), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e30 * .e43 * .e5 * .e29 * .e57 * .e12/((.e33 * .e12)^2 *
    (1 + .e28 + .e29) * (.e8) * .e11)))), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e30 * .e43 * .e6 * .e29 * .e60 * .e12/((.e33 * .e12)^2 *
    (1 + .e28 + .e29) * (.e8) * .e11)))), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    nZHvar)] <- crossprod(Xvar * wHvar * (-(S * (2 * ((2 * (.e16 * .e23/.e10) -
    .e33)/(.e33 * (1 + .e28 + .e29))^2) + 2/(.e33 * (1 + .e28 + .e29)^2)) * .e43 *
    .e28 * .e29/((.e8) * .e11))), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * S * .e69 * .e43 *
    .e29/((.e8) * .e11), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e3 * (S * (0.5 * (S * .e17 * (S * (0.5 * (S * .e25 * epsilon2/(.e8)^2) -
    .e47 * .e24) * epsilon2 - 2 * (.e25/(.e8))) * epsilon2/(.e8)^2) - (0.5 *
    (S^2 * .e47 * .e17 * epsilon2^2/(.e8)^2) - (((0.5 * (.e3/(.e8)) + 1 - 0.5 *
    (0.5 * (1 - .e3/(.e8)) + .e3/(.e8))) * (1 - .e3/(.e8)) * .e4/.e15 + (2 -
    2 * (.e46^2 * .e3 * (.e8)/((.e8) * .e15)^2)) * .e15)/((.e8) * .e15)^2 + S^2 *
    .e47^2 * .e3 * epsilon2^2/((.e8) * .e15)) * .e17) * .e24) * epsilon2 - 0.5 *
    ((S * (0.5 * .e45 - .e47 * .e24 * .e17) * epsilon2 - .e17 * .e25/(.e8))/(.e8))) +
    S * (0.5 * .e45 - .e47 * .e24 * .e17) * epsilon2 - 0.5 * (.e17 * .e25/(.e8)))/.e49 -
    (0.5 * (.e33 * (1 + .e28 + .e29)/.e11) + 2 * (.e29 * .e48)) * .e3 * .e48/.e49^2) *
    .e3 * .e29), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((S * (((((0.5 * ((1 - .e3/(.e8)) * .e4) - S^2 * .e50 * .e47 *
    .e3 * epsilon2^2)/(.e8) + 0.5 * ((.e3/(.e8) - 1) * .e4/(.e8) + 1 - 0.5 *
    ((1 - .e3/(.e8)) * (1 - .e4/(.e8))))) * .e17/.e15 + 0.5 * (S^2 * .e50 * .e17 *
    epsilon2^2/(.e8)^2)) * .e3 + .e50 * (1 - 2 * (.e46 * .e3 * (.e8) * .e15/((.e8) *
    .e15)^2)) * .e17) * .e24/((.e8) * .e15)^2 + 0.5 * (S * .e17 * (S * (0.5 *
    (S * .e25 * epsilon2/(.e8)^2) - .e47 * .e24) * epsilon2 - 2 * (.e25/(.e8))) *
    epsilon2/(.e8)^2)) * epsilon2 - 0.5 * ((S * (0.5 * .e45 - .e47 * .e24 * .e17) *
    epsilon2 - .e17 * .e25/(.e8))/(.e8)))/.e49 - (0.5 * (.e33 * (1 + .e28 + .e29)/.e11) +
    2 * (.e29 * .e48)) * .e52/.e49^2) * .e3 * .e4 * .e29), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e30 * .e53 * .e3 *
    (.e9) * .e29 * .e48 * .e12/((.e33 * (.e9) * .e12)^2 * (1 + .e28 + .e29) *
    .e11)))), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e30 * .e3 * .e5 * .e29 *
    .e48 * .e57 * .e12/((.e33 * .e12)^2 * (1 + .e28 + .e29) * .e11)))), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e30 * .e3 * .e6 * .e29 *
    .e60 * .e48 * .e12/((.e33 * .e12)^2 * (1 + .e28 + .e29) * .e11)))), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-((2 * ((2 * (.e16 *
    .e23/.e10) - .e33)/(.e33 * (1 + .e28 + .e29))^2) + 2/(.e33 * (1 + .e28 +
    .e29)^2)) * .e3 * .e28 * .e29 * .e48/.e11)), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar + 3 * nuZUvar +
      3 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * .e69 * .e3 *
    .e29 * .e48/.e11, Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 * (.e4/(.e8)) -
    0.5 * (0.5 * (1 - .e4/(.e8)) + .e4/(.e8))) * (1 - .e4/(.e8)) + S^2 * .e50^2 *
    .e3 * .e4 * epsilon2^2/(((.e8) * .e15)^2 * (.e8))) * .e17 * .e3/.e15 + ((0.5 *
    (S^2 * .e17 * epsilon2^2/(.e8)^2) - 2 * (.e50 * .e17 * (.e8) * .e15/((.e8) *
    .e15)^2)) * .e4 + .e17) * .e50) * .e24 * .e3/((.e8) * .e15)^2 + S * (0.5 *
    (.e4 * (S * (.e50 * .e24 * .e3/((.e8) * .e15)^2 + 0.5 * (S * .e25 * epsilon2/(.e8)^2)) *
      epsilon2 - 2 * (.e25/(.e8)))) + 0.5 * .e25) * .e17 * epsilon2/(.e8)^2) *
    epsilon2 - (0.5 * (.e17 * .e25) + 0.5 * (.e4 * (S * .e51 * epsilon2 - .e17 *
    .e25/(.e8))))/(.e8))/.e49 - (0.5 * (.e33 * (1 + .e28 + .e29)/.e11) + 2 *
    (.e29 * .e52)) * .e4 * .e52/.e49^2) * .e4 * .e29), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e30 * .e53 * (.e9) *
    .e4 * .e29 * .e52 * .e12/((.e33 * (.e9) * .e12)^2 * (1 + .e28 + .e29) * .e11)))),
    Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e30 * .e5 * .e4 * .e29 *
    .e52 * .e57 * .e12/((.e33 * .e12)^2 * (1 + .e28 + .e29) * .e11)))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e30 * .e4 * .e6 * .e29 *
    .e52 * .e60 * .e12/((.e33 * .e12)^2 * (1 + .e28 + .e29) * .e11)))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-((2 * ((2 * (.e16 *
    .e23/.e10) - .e33)/(.e33 * (1 + .e28 + .e29))^2) + 2/(.e33 * (1 + .e28 +
    .e29)^2)) * .e4 * .e28 * .e29 * .e52/.e11)), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * .e69 *
    .e4 * .e29 * .e52/.e11, Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (((.e18 * (S * (.e26 *
    .e5/.e14 + S * .e27 * epsilon3) * epsilon3/(.e9) - .e27) + S * .e26 * (.e18 *
    .e5/.e6 + .e18) * .e5 * epsilon3/((.e9) * .e14))/(.e33 * (.e9) * .e12) -
    2 * (.e30 * .e53^2/(.e33 * (.e9) * .e12)^2)) * .e30), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * ((.e56 * .e26 *
    .e18 + (S * ((0.5 * (S * (.e26 * .e5/.e14 + S * .e27 * epsilon3) * epsilon3/(.e9) -
    .e27) - 0.5 * .e27) * .e18/(.e9) - S * .e56 * .e26 * (.e18 * .e5/.e6 + .e18) *
    epsilon3) * epsilon3 - 0.5 * (.e53/(.e9)))/(.e9))/(.e33 * .e12) - 2 * (.e30 *
    .e53 * .e57/((.e33 * .e12)^2 * (.e9)))) * .e30 * .e5), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * (((S * ((0.5 *
    (S * (.e26 * .e5/.e14 + S * .e27 * epsilon3) * epsilon3/(.e9) - .e27) - 0.5 *
    .e27) * .e18/(.e9) + S * .e58 * .e26 * (.e18 * .e5/.e6 + .e18) * .e5 * epsilon3/((.e9) *
    .e14)^2) * epsilon3 - 0.5 * (.e53/(.e9)))/(.e9) - .e58 * .e26 * .e18 * .e5/((.e9) *
    .e14)^2)/(.e33 * .e12) - 2 * (.e30 * .e53 * .e60/((.e33 * .e12)^2 * (.e9)))) *
    .e30 * .e6), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * (-(S * .e30 *
    .e71 * .e53 * .e28/((.e9) * .e12))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e30 * .e73 * .e53 * .e29/((.e9) * .e12))), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((.e5 * (S * (0.5 *
    (S * .e18 * (S * (0.5 * (S * .e27 * epsilon3/(.e9)^2) - .e56 * .e26) * epsilon3 -
      2 * (.e27/(.e9))) * epsilon3/(.e9)^2) - (0.5 * (S^2 * .e56 * .e18 * epsilon3^2/(.e9)^2) -
    (((0.5 * (.e5/(.e9)) + 1 - 0.5 * (0.5 * (1 - .e5/(.e9)) + .e5/(.e9))) * (1 -
      .e5/(.e9)) * .e6/.e14 + (2 - 2 * (.e55^2 * .e5 * (.e9)/((.e9) * .e14)^2)) *
      .e14)/((.e9) * .e14)^2 + S^2 * .e56^2 * .e5 * epsilon3^2/((.e9) * .e14)) *
      .e18) * .e26) * epsilon3 - 0.5 * ((S * (0.5 * .e54 - .e56 * .e26 * .e18) *
    epsilon3 - .e18 * .e27/(.e9))/(.e9))) + S * (0.5 * .e54 - .e56 * .e26 * .e18) *
    epsilon3 - 0.5 * (.e18 * .e27/(.e9)))/(.e33 * .e12) - (0.5 * (.e33/.e12) +
    2 * (.e30 * .e57)) * .e5 * .e57/(.e33 * .e12)^2) * .e30 * .e5), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 *
    ((1 - .e5/(.e9)) * .e6) - S^2 * .e58 * .e56 * .e5 * epsilon3^2)/(.e9) + 0.5 *
    ((.e5/(.e9) - 1) * .e6/(.e9) + 1 - 0.5 * ((1 - .e5/(.e9)) * (1 - .e6/(.e9))))) *
    .e18/.e14 + 0.5 * (S^2 * .e58 * .e18 * epsilon3^2/(.e9)^2)) * .e5 + .e58 *
    (1 - 2 * (.e55 * .e5 * (.e9) * .e14/((.e9) * .e14)^2)) * .e18) * .e26/((.e9) *
    .e14)^2 + 0.5 * (S * .e18 * (S * (0.5 * (S * .e27 * epsilon3/(.e9)^2) - .e56 *
    .e26) * epsilon3 - 2 * (.e27/(.e9))) * epsilon3/(.e9)^2)) * epsilon3 - 0.5 *
    ((S * (0.5 * .e54 - .e56 * .e26 * .e18) * epsilon3 - .e18 * .e27/(.e9))/(.e9)))/(.e33 *
    .e12) - (0.5 * (.e33/.e12) + 2 * (.e30 * .e57)) * .e60/(.e33 * .e12)^2) *
    .e30 * .e5 * .e6), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-(.e30 * .e71 *
    .e5 * .e28 * .e57/.e12)), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e30 *
    .e73 * .e5 * .e29 * .e57/.e12)), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e6/(.e9)) - 0.5 * (0.5 * (1 - .e6/(.e9)) + .e6/(.e9))) * (1 - .e6/(.e9)) +
    S^2 * .e58^2 * .e5 * .e6 * epsilon3^2/(((.e9) * .e14)^2 * (.e9))) * .e18 *
    .e5/.e14 + ((0.5 * (S^2 * .e18 * epsilon3^2/(.e9)^2) - 2 * (.e58 * .e18 *
    (.e9) * .e14/((.e9) * .e14)^2)) * .e6 + .e18) * .e58) * .e26 * .e5/((.e9) *
    .e14)^2 + S * (0.5 * (.e6 * (S * (.e58 * .e26 * .e5/((.e9) * .e14)^2 + 0.5 *
    (S * .e27 * epsilon3/(.e9)^2)) * epsilon3 - 2 * (.e27/(.e9)))) + 0.5 * .e27) *
    .e18 * epsilon3/(.e9)^2) * epsilon3 - (0.5 * (.e18 * .e27) + 0.5 * (.e6 *
    (S * .e59 * epsilon3 - .e18 * .e27/(.e9))))/(.e9))/(.e33 * .e12) - (0.5 *
    (.e33/.e12) + 2 * (.e30 * .e60)) * .e6 * .e60/(.e33 * .e12)^2) * .e30 * .e6),
    vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-(.e30 * .e71 *
    .e6 * .e28 * .e60/.e12)), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e30 *
    .e73 * .e6 * .e29 * .e60/.e12)), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + nZHvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + nZHvar)] <- crossprod(Zvar * wHvar * ((1 - .e28/(1 +
    .e28 + .e29))/(.e33 * (1 + .e28 + .e29)) - 2 * (.e16 * .e28 * .e23/((.e33 *
    (1 + .e28 + .e29))^2 * .e10))) * (2 * (.e16 * .e23/.e10) - .e33) * .e28,
    Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + nZHvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-(((2 * (.e16 * .e23/.e10) - .e33)/(.e33 * (1 + .e28 + .e29)^2) +
    2 * ((2 * (.e17 * .e25/.e11) - .e33) * .e16 * .e23/((.e33 * (1 + .e28 + .e29))^2 *
      .e10))) * .e28 * .e29)), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar + 2 * nZHvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    nZHvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar)] <- crossprod(Zvar *
    wHvar * ((1 - .e29/(1 + .e28 + .e29))/(.e33 * (1 + .e28 + .e29)) - 2 * (.e17 *
    .e29 * .e25/((.e33 * (1 + .e28 + .e29))^2 * .e11))) * (2 * (.e17 * .e25/.e11) -
    .e33) * .e29, Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for lcm 3 classes halfnormal-normal distribution
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
LCM3ChnormAlgOpt <- function(start, olsParam, dataTable, S, wHvar,
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csLCMfhalfnorm3C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(cLCMhalfnormlike3C(startVal, nXvar = nXvar,
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
  cat("LCM 3 Classes Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cLCMhalfnormlike3C(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike3C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cLCMhalfnormlike3C,
    grad = cgradLCMhalfnormlike3C, hess = chessLCMhalfnormlike3C,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cLCMhalfnormlike3C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessLCMhalfnormlike3C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradLCMhalfnormlike3C(mleObj$par,
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
      mleObj$hessian <- chessLCMhalfnormlike3C(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessLCMhalfnormlike3C(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cLCMhalfnormlike3C(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradLCMhalfnormlike3C(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) initHalf = initHalf))
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for lcmcross 3 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @param level level for confidence interval
#' @noRd
cLCM3Chalfnormeff <- function(object, level) {
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
  theta1 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar + 1):(3 * object$nXvar +
    3 * object$nuZUvar + 3 * object$nvZVvar + 2 * object$nZHvar)]
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
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
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
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
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
    teBC_c3 <- exp(-mustar3 + 1/2 * sigmastar3^2) * pnorm(mustar3/sigmastar3 -
      sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, ifelse(Group_c ==
      2, teBC_c2, teBC_c3))
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    effBC_c3 <- ifelse(Group_c == 3, teBC_c3, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) *
      pnorm(mustar1/sigmastar1 + sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) *
      pnorm(mustar2/sigmastar2 + sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c3 <- exp(mustar3 + 1/2 * sigmastar3^2) *
      pnorm(mustar3/sigmastar3 + sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      ifelse(Group_c == 2, teBC_reciprocal_c2, teBC_reciprocal_c3))
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    ReffBC_c3 <- ifelse(Group_c == 3, teBC_reciprocal_c3,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
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
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
      u_c3 = u_c3, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2,
      ineff_c3 = ineff_c3)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal effects for for lcmcross 3 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @noRd
cmargLCM3Chalfnorm_Eu <- function(object) {
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
  theta1 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar + 1):(3 * object$nXvar +
    3 * object$nuZUvar + 3 * object$nvZVvar + 2 * object$nZHvar)]
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
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
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
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu3/2) * dnorm(0), ncol = 1))
  colnames(margEff_c3) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c3")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      ifelse(Group_c == 2, margEff_c2[, c], margEff_c3[,
        c]))
  }
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2,
    margEff_c3)
  return(margEff)
}

cmargLCM3Chalfnorm_Vu <- function(object) {
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
  theta1 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
    3 * object$nvZVvar + object$nZHvar + 1):(3 * object$nXvar +
    3 * object$nuZUvar + 3 * object$nvZVvar + 2 * object$nZHvar)]
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
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
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
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu3) * (1 - (dnorm(0)/pnorm(0))^2),
    ncol = 1))
  colnames(margEff_c3) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c3")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      ifelse(Group_c == 2, margEff_c2[, c], margEff_c3[,
        c]))
  }
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2,
    margEff_c3)
  return(margEff)
}
