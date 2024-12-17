################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Latent Class Stochastic Frontier Analysis                             #
# Number of Classes: 5L                                                        #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for lcm 5 classes halfnormal-normal distribution
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
cLCMhalfnormlike5C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  beta5 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  delta5 <- parm[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 4 * nvZVvar)]
  phi5 <- parm[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar)]
  theta1 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar)]
  theta2 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 *
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    3 * nZHvar)]
  theta4 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 *
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    4 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon5 <- Yvar - as.numeric(crossprod(matrix(beta5), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  mustar5 <- -exp(Wu5) * S * epsilon5/(exp(Wu5) + exp(Wv5))
  sigmastar5 <- sqrt(exp(Wu5) * exp(Wv5)/(exp(Wu5) + exp(Wv5)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Pi5 <- 2/sqrt(exp(Wu5) + exp(Wv5)) * dnorm(S * epsilon5/sqrt(exp(Wu5) +
    exp(Wv5))) * pnorm(mustar5/sigmastar5)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc4 <- exp(Wz4)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc5 <- 1 - Probc1 - Probc2 - Probc3 - Probc4
  ll <- log(Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 + Probc4 * Pi4 +
    Probc5 * Pi5)
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for lcm 5 classes halfnormal-normal distribution
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
csLCMfhalfnorm5C <- function(olsObj, epsiRes, nXvar, nuZUvar,
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
  StartVal <- c(Esti[1:(nXvar + 1)], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0,
    nvZVvar - 1), 0.98 * Esti[1:nXvar], Esti[nXvar + 1],
    if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar + 2],
    if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:nXvar],
    Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1),
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1),
    0.98 * Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
      1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
      1) rep(0, nvZVvar - 1), 0.98 * Esti[1:nXvar], Esti[nXvar +
      1], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar +
      2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0,
      4 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)),
    paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)),
    paste0("Zv_", colnames(vHvar)), paste0("Cl1_", colnames(Zvar)),
    paste0("Cl2_", colnames(Zvar)), paste0("Cl3_", colnames(Zvar)),
    paste0("Cl4_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for lcm 5 classes halfnormal-normal distribution
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
cgradLCMhalfnormlike5C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  beta5 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  delta5 <- parm[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 4 * nvZVvar)]
  phi5 <- parm[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar)]
  theta1 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar)]
  theta2 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 *
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    3 * nZHvar)]
  theta4 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 *
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    4 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon5 <- Yvar - as.numeric(crossprod(matrix(beta5), t(Xvar)))
  .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wu3)
  .e6 <- exp(Wv3)
  .e7 <- exp(Wu4)
  .e8 <- exp(Wv4)
  .e9 <- exp(Wu5)
  .e10 <- exp(Wv5)
  .e11 <- exp(Wz1)
  .e12 <- exp(Wz2)
  .e13 <- exp(Wz3)
  .e14 <- exp(Wz4)
  .e15 <- .e1 + .e2
  .e16 <- .e3 + .e4
  .e17 <- .e5 + .e6
  .e18 <- .e7 + .e8
  .e19 <- .e9 + .e10
  .e20 <- sqrt(.e1 * .e2/(.e15))
  .e21 <- sqrt(.e3 * .e4/(.e16))
  .e22 <- sqrt(.e5 * .e6/(.e17))
  .e23 <- sqrt(.e7 * .e8/(.e18))
  .e24 <- sqrt(.e9 * .e10/(.e19))
  .e25 <- (S * .e1 * epsilon1/((.e15) * .e20))
  .e26 <- (S * .e3 * epsilon2/((.e16) * .e21))
  .e27 <- (S * .e5 * epsilon3/((.e17) * .e22))
  .e28 <- (S * .e7 * epsilon4/((.e18) * .e23))
  .e29 <- (S * .e9 * epsilon5/((.e19) * .e24))
  .e30 <- pnorm(-.e25)
  .e31 <- pnorm(-.e26)
  .e32 <- pnorm(-.e27)
  .e33 <- pnorm(-.e28)
  .e34 <- pnorm(-.e29)
  .e35 <- dnorm(S * epsilon1/sqrt(.e15))
  .e36 <- dnorm(S * epsilon2/sqrt(.e16))
  .e37 <- dnorm(S * epsilon3/sqrt(.e17))
  .e38 <- dnorm(S * epsilon4/sqrt(.e18))
  .e39 <- dnorm(S * epsilon5/sqrt(.e19))
  .e40 <- dnorm(-.e25)
  .e41 <- dnorm(-.e26)
  .e42 <- dnorm(-.e27)
  .e43 <- dnorm(-.e28)
  .e44 <- dnorm(-.e29)
  .e45 <- (.e35 * .e11 * .e30/sqrt(.e15))
  .e46 <- (.e36 * .e12 * .e31/sqrt(.e16))
  .e47 <- (.e37 * .e13 * .e32/sqrt(.e17))
  .e48 <- (.e38 * .e14 * .e33/sqrt(.e18))
  .e49 <- (1 - (.e11 + .e12 + .e13 + .e14)/(1 + .e11 + .e12 + .e13 + .e14))
  .e50 <- (.e49 * .e39 * .e34/sqrt(.e19))
  .e51 <- (.e40 * .e35 * .e1/.e20 + S * .e35 * .e30 * epsilon1)
  .e52 <- (2 * .e45 + 2 * .e46 + 2 * .e47 + 2 * .e48)
  .e53 <- (.e52/(1 + .e11 + .e12 + .e13 + .e14) + 2 * .e50)
  .e54 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e15) * sqrt(.e15))
  .e55 <- (0.5 * ((1 - .e1/(.e15)) * .e2/.e20) + .e20)
  .e56 <- (1/((.e15) * .e20) - .e55 * .e1/((.e15) * .e20)^2)
  .e57 <- (0.5 * (S * .e35 * .e30 * epsilon1/(.e15)^2) - .e56 * .e40 * .e35)
  .e58 <- (S * .e57 * epsilon1 - 0.5 * (.e35 * .e30/(.e15)))
  .e59 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))
  .e60 <- (S * .e35 * .e30 * epsilon1/(.e15)^2)
  .e61 <- (0.5 * ((1 - .e2/(.e15)) * .e1/.e20) + .e20)
  .e62 <- (.e61 * .e40 * .e35 * .e1/((.e15) * .e20)^2 + 0.5 * .e60)
  .e63 <- (S * .e62 * epsilon1 - 0.5 * (.e35 * .e30/(.e15)))
  .e64 <- (.e41 * .e36 * .e3/.e21 + S * .e36 * .e31 * epsilon2)
  .e65 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e16) * sqrt(.e16))
  .e66 <- (0.5 * ((1 - .e3/(.e16)) * .e4/.e21) + .e21)
  .e67 <- (1/((.e16) * .e21) - .e66 * .e3/((.e16) * .e21)^2)
  .e68 <- (0.5 * (S * .e36 * .e31 * epsilon2/(.e16)^2) - .e67 * .e41 * .e36)
  .e69 <- (S * .e68 * epsilon2 - 0.5 * (.e36 * .e31/(.e16)))
  .e70 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e16))
  .e71 <- (S * .e36 * .e31 * epsilon2/(.e16)^2)
  .e72 <- (0.5 * ((1 - .e4/(.e16)) * .e3/.e21) + .e21)
  .e73 <- (.e72 * .e41 * .e36 * .e3/((.e16) * .e21)^2 + 0.5 * .e71)
  .e74 <- (S * .e73 * epsilon2 - 0.5 * (.e36 * .e31/(.e16)))
  .e75 <- (.e42 * .e37 * .e5/.e22 + S * .e37 * .e32 * epsilon3)
  .e76 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e17) * sqrt(.e17))
  .e77 <- (0.5 * ((1 - .e5/(.e17)) * .e6/.e22) + .e22)
  .e78 <- (1/((.e17) * .e22) - .e77 * .e5/((.e17) * .e22)^2)
  .e79 <- (0.5 * (S * .e37 * .e32 * epsilon3/(.e17)^2) - .e78 * .e42 * .e37)
  .e80 <- (S * .e79 * epsilon3 - 0.5 * (.e37 * .e32/(.e17)))
  .e81 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e17))
  .e82 <- (S * .e37 * .e32 * epsilon3/(.e17)^2)
  .e83 <- (0.5 * ((1 - .e6/(.e17)) * .e5/.e22) + .e22)
  .e84 <- (.e83 * .e42 * .e37 * .e5/((.e17) * .e22)^2 + 0.5 * .e82)
  .e85 <- (S * .e84 * epsilon3 - 0.5 * (.e37 * .e32/(.e17)))
  .e86 <- (.e43 * .e38 * .e7/.e23 + S * .e38 * .e33 * epsilon4)
  .e87 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e18) * sqrt(.e18))
  .e88 <- (0.5 * ((1 - .e7/(.e18)) * .e8/.e23) + .e23)
  .e89 <- (1/((.e18) * .e23) - .e88 * .e7/((.e18) * .e23)^2)
  .e90 <- (0.5 * (S * .e38 * .e33 * epsilon4/(.e18)^2) - .e89 * .e43 * .e38)
  .e91 <- (S * .e90 * epsilon4 - 0.5 * (.e38 * .e33/(.e18)))
  .e92 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e18))
  .e93 <- (S * .e38 * .e33 * epsilon4/(.e18)^2)
  .e94 <- (0.5 * ((1 - .e8/(.e18)) * .e7/.e23) + .e23)
  .e95 <- (.e94 * .e43 * .e38 * .e7/((.e18) * .e23)^2 + 0.5 * .e93)
  .e96 <- (S * .e95 * epsilon4 - 0.5 * (.e38 * .e33/(.e18)))
  .e97 <- (.e44 * .e39 * .e9/.e24 + S * .e39 * .e34 * epsilon5)
  .e98 <- (.e53 * (.e19) * sqrt(.e19))
  .e99 <- (0.5 * ((1 - .e9/(.e19)) * .e10/.e24) + .e24)
  .e100 <- (1/((.e19) * .e24) - .e99 * .e9/((.e19) * .e24)^2)
  .e101 <- (S * .e39 * .e34 * epsilon5/(.e19)^2)
  .e102 <- (S * (0.5 * .e101 - .e100 * .e44 * .e39) * epsilon5 - 0.5 * (.e39 *
    .e34/(.e19)))
  .e103 <- (0.5 * ((1 - .e10/(.e19)) * .e9/.e24) + .e24)
  .e104 <- (.e103 * .e44 * .e39 * .e9/((.e19) * .e24)^2 + 0.5 * .e101)
  .e105 <- (S * .e104 * epsilon5 - 0.5 * (.e39 * .e34/(.e19)))
  .e106 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14))
  .e107 <- (2 * (.e35 * .e30/sqrt(.e15)) - .e53)
  .e108 <- (2 * (.e36 * .e31/sqrt(.e16)) - .e53)
  .e109 <- (2 * (.e37 * .e32/sqrt(.e17)) - .e53)
  .e110 <- (2 * (.e38 * .e33/sqrt(.e18)) - .e53)
  gradll <- cbind(2 * (S * Xvar * .e51 * .e11/.e54), 2 * (uHvar * .e1 * .e11 *
    .e58/.e59), 2 * (vHvar * .e2 * .e11 * .e63/.e59), 2 * (S * Xvar * .e64 *
    .e12/.e65), 2 * (uHvar * .e3 * .e12 * .e69/.e70), 2 * (vHvar * .e4 * .e12 *
    .e74/.e70), 2 * (S * Xvar * .e75 * .e13/.e76), 2 * (uHvar * .e5 * .e13 *
    .e80/.e81), 2 * (vHvar * .e6 * .e13 * .e85/.e81), 2 * (S * Xvar * .e86 *
    .e14/.e87), 2 * (uHvar * .e7 * .e14 * .e91/.e92), 2 * (vHvar * .e8 * .e14 *
    .e96/.e92), 2 * (S * Xvar * .e49 * .e97/.e98), 2 * (uHvar * .e49 * .e9 *
    .e102/(.e53 * sqrt(.e19))), 2 * (vHvar * .e49 * .e10 * .e105/(.e53 * sqrt(.e19))),
    Zvar * .e107 * .e11/.e106, Zvar * .e108 * .e12/.e106, Zvar * .e109 * .e13/.e106,
    Zvar * .e110 * .e14/.e106)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for lcm 5 classes halfnormal-normal distribution
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
chessLCMhalfnormlike5C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  beta5 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  delta5 <- parm[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 4 * nvZVvar)]
  phi5 <- parm[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar)]
  theta1 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar)]
  theta2 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 *
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    3 * nZHvar)]
  theta4 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 *
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    4 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon5 <- Yvar - as.numeric(crossprod(matrix(beta5), t(Xvar)))
  .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wu3)
  .e6 <- exp(Wv3)
  .e7 <- exp(Wu4)
  .e8 <- exp(Wv4)
  .e9 <- exp(Wu5)
  .e10 <- exp(Wv5)
  .e11 <- exp(Wz1)
  .e12 <- exp(Wz2)
  .e13 <- exp(Wz3)
  .e14 <- exp(Wz4)
  .e15 <- .e1 + .e2
  .e16 <- .e3 + .e4
  .e17 <- .e5 + .e6
  .e18 <- .e7 + .e8
  .e19 <- .e9 + .e10
  .e20 <- sqrt(.e1 * .e2/(.e15))
  .e21 <- sqrt(.e3 * .e4/(.e16))
  .e22 <- sqrt(.e5 * .e6/(.e17))
  .e23 <- sqrt(.e7 * .e8/(.e18))
  .e24 <- sqrt(.e9 * .e10/(.e19))
  .e25 <- (S * .e1 * epsilon1/((.e15) * .e20))
  .e26 <- (S * .e3 * epsilon2/((.e16) * .e21))
  .e27 <- (S * .e5 * epsilon3/((.e17) * .e22))
  .e28 <- (S * .e7 * epsilon4/((.e18) * .e23))
  .e29 <- (S * .e9 * epsilon5/((.e19) * .e24))
  .e30 <- pnorm(-.e25)
  .e31 <- pnorm(-.e26)
  .e32 <- pnorm(-.e27)
  .e33 <- pnorm(-.e28)
  .e34 <- pnorm(-.e29)
  .e35 <- dnorm(S * epsilon1/sqrt(.e15))
  .e36 <- dnorm(S * epsilon2/sqrt(.e16))
  .e37 <- dnorm(S * epsilon3/sqrt(.e17))
  .e38 <- dnorm(S * epsilon4/sqrt(.e18))
  .e39 <- dnorm(S * epsilon5/sqrt(.e19))
  .e40 <- dnorm(-.e25)
  .e41 <- dnorm(-.e26)
  .e42 <- dnorm(-.e27)
  .e43 <- dnorm(-.e28)
  .e44 <- dnorm(-.e29)
  .e45 <- (.e35 * .e11 * .e30/sqrt(.e15))
  .e46 <- (.e36 * .e12 * .e31/sqrt(.e16))
  .e47 <- (.e37 * .e13 * .e32/sqrt(.e17))
  .e48 <- (.e38 * .e14 * .e33/sqrt(.e18))
  .e49 <- (1 - (.e11 + .e12 + .e13 + .e14)/(1 + .e11 + .e12 + .e13 + .e14))
  .e50 <- (.e49 * .e39 * .e34/sqrt(.e19))
  .e51 <- (.e40 * .e35 * .e1/.e20 + S * .e35 * .e30 * epsilon1)
  .e52 <- (2 * .e45 + 2 * .e46 + 2 * .e47 + 2 * .e48)
  .e53 <- (.e52/(1 + .e11 + .e12 + .e13 + .e14) + 2 * .e50)
  .e54 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e15) * sqrt(.e15))
  .e55 <- (0.5 * ((1 - .e1/(.e15)) * .e2/.e20) + .e20)
  .e56 <- (1/((.e15) * .e20) - .e55 * .e1/((.e15) * .e20)^2)
  .e57 <- (0.5 * (S * .e35 * .e30 * epsilon1/(.e15)^2) - .e56 * .e40 * .e35)
  .e58 <- (S * .e57 * epsilon1 - 0.5 * (.e35 * .e30/(.e15)))
  .e59 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))
  .e60 <- (S * .e35 * .e30 * epsilon1/(.e15)^2)
  .e61 <- (0.5 * ((1 - .e2/(.e15)) * .e1/.e20) + .e20)
  .e62 <- (.e61 * .e40 * .e35 * .e1/((.e15) * .e20)^2 + 0.5 * .e60)
  .e63 <- (S * .e62 * epsilon1 - 0.5 * (.e35 * .e30/(.e15)))
  .e64 <- (.e41 * .e36 * .e3/.e21 + S * .e36 * .e31 * epsilon2)
  .e65 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e16) * sqrt(.e16))
  .e66 <- (0.5 * ((1 - .e3/(.e16)) * .e4/.e21) + .e21)
  .e67 <- (1/((.e16) * .e21) - .e66 * .e3/((.e16) * .e21)^2)
  .e68 <- (0.5 * (S * .e36 * .e31 * epsilon2/(.e16)^2) - .e67 * .e41 * .e36)
  .e69 <- (S * .e68 * epsilon2 - 0.5 * (.e36 * .e31/(.e16)))
  .e70 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e16))
  .e71 <- (S * .e36 * .e31 * epsilon2/(.e16)^2)
  .e72 <- (0.5 * ((1 - .e4/(.e16)) * .e3/.e21) + .e21)
  .e73 <- (.e72 * .e41 * .e36 * .e3/((.e16) * .e21)^2 + 0.5 * .e71)
  .e74 <- (S * .e73 * epsilon2 - 0.5 * (.e36 * .e31/(.e16)))
  .e75 <- (.e42 * .e37 * .e5/.e22 + S * .e37 * .e32 * epsilon3)
  .e76 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e17) * sqrt(.e17))
  .e77 <- (0.5 * ((1 - .e5/(.e17)) * .e6/.e22) + .e22)
  .e78 <- (1/((.e17) * .e22) - .e77 * .e5/((.e17) * .e22)^2)
  .e79 <- (0.5 * (S * .e37 * .e32 * epsilon3/(.e17)^2) - .e78 * .e42 * .e37)
  .e80 <- (S * .e79 * epsilon3 - 0.5 * (.e37 * .e32/(.e17)))
  .e81 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e17))
  .e82 <- (S * .e37 * .e32 * epsilon3/(.e17)^2)
  .e83 <- (0.5 * ((1 - .e6/(.e17)) * .e5/.e22) + .e22)
  .e84 <- (.e83 * .e42 * .e37 * .e5/((.e17) * .e22)^2 + 0.5 * .e82)
  .e85 <- (S * .e84 * epsilon3 - 0.5 * (.e37 * .e32/(.e17)))
  .e86 <- (.e43 * .e38 * .e7/.e23 + S * .e38 * .e33 * epsilon4)
  .e87 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * (.e18) * sqrt(.e18))
  .e88 <- (0.5 * ((1 - .e7/(.e18)) * .e8/.e23) + .e23)
  .e89 <- (1/((.e18) * .e23) - .e88 * .e7/((.e18) * .e23)^2)
  .e90 <- (0.5 * (S * .e38 * .e33 * epsilon4/(.e18)^2) - .e89 * .e43 * .e38)
  .e91 <- (S * .e90 * epsilon4 - 0.5 * (.e38 * .e33/(.e18)))
  .e92 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e18))
  .e93 <- (S * .e38 * .e33 * epsilon4/(.e18)^2)
  .e94 <- (0.5 * ((1 - .e8/(.e18)) * .e7/.e23) + .e23)
  .e95 <- (.e94 * .e43 * .e38 * .e7/((.e18) * .e23)^2 + 0.5 * .e93)
  .e96 <- (S * .e95 * epsilon4 - 0.5 * (.e38 * .e33/(.e18)))
  .e97 <- (.e44 * .e39 * .e9/.e24 + S * .e39 * .e34 * epsilon5)
  .e98 <- (.e53 * (.e19) * sqrt(.e19))
  .e99 <- (0.5 * ((1 - .e9/(.e19)) * .e10/.e24) + .e24)
  .e100 <- (1/((.e19) * .e24) - .e99 * .e9/((.e19) * .e24)^2)
  .e101 <- (S * .e39 * .e34 * epsilon5/(.e19)^2)
  .e102 <- (S * (0.5 * .e101 - .e100 * .e44 * .e39) * epsilon5 - 0.5 * (.e39 *
    .e34/(.e19)))
  .e103 <- (0.5 * ((1 - .e10/(.e19)) * .e9/.e24) + .e24)
  .e104 <- (.e103 * .e44 * .e39 * .e9/((.e19) * .e24)^2 + 0.5 * .e101)
  .e105 <- (S * .e104 * epsilon5 - 0.5 * (.e39 * .e34/(.e19)))
  .e106 <- (.e53 * (1 + .e11 + .e12 + .e13 + .e14))
  .e107 <- (2 * (.e35 * .e30/sqrt(.e15)) - .e53)
  .e108 <- (2 * (.e36 * .e31/sqrt(.e16)) - .e53)
  .e109 <- (2 * (.e37 * .e32/sqrt(.e17)) - .e53)
  .e110 <- (2 * (.e38 * .e33/sqrt(.e18)) - .e53)
  .e111 <- (S * (.e40 * .e1/.e20 + S * .e30 * epsilon1) * epsilon1/(.e15) - .e30)
  .e112 <- (2 * (.e108/.e106^2) + 2/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2))
  .e113 <- (2 * (.e109/.e106^2) + 2/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2))
  .e114 <- (2 * (.e110/.e106^2) + 2/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2))
  .e115 <- 2 * (.e107/.e106^2) + 2/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2)
  .e116 <- ((2 - 2 * (.e12/(1 + .e11 + .e12 + .e13 + .e14)))/.e106 - 2 * (.e108 *
    .e12/.e106^2))
  .e117 <- ((2 - 2 * (.e13/(1 + .e11 + .e12 + .e13 + .e14)))/.e106 - 2 * (.e109 *
    .e13/.e106^2))
  .e118 <- ((2 - 2 * (.e14/(1 + .e11 + .e12 + .e13 + .e14)))/.e106 - 2 * (.e110 *
    .e14/.e106^2))
  .e119 <- (2 * ((1 + .e11 + .e12 + .e13 + .e14) * .e107/.e106^2) + 2/.e106)
  .e120 <- (2 * ((1 + .e11 + .e12 + .e13 + .e14) * .e108/.e106^2) + 2/.e106)
  .e121 <- (2 * ((1 + .e11 + .e12 + .e13 + .e14) * .e110/.e106^2) + 2/.e106)
  .e122 <- (2 * ((1 + .e11 + .e12 + .e13 + .e14) * .e109/.e106^2) + 2/.e106)
  hessll <- matrix(nrow = (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar),
    ncol = (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar))
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * wHvar * 2 * (((.e35 * .e111 + S *
    .e40 * (.e35 * .e1/.e2 + .e35) * .e1 * epsilon1/((.e15) * .e20))/.e54 - 2 *
    (.e51^2 * .e11/.e54^2)) * .e11), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(Xvar * wHvar * 2 *
    (S * ((.e56 * .e40 * .e35 + (S * ((0.5 * .e111 - 0.5 * .e30) * .e35/(.e15) -
      S * .e56 * .e40 * (.e35 * .e1/.e2 + .e35) * epsilon1) * epsilon1 - 0.5 *
      (.e51/(.e15)))/(.e15))/.e59 - 2 * (.e51 * .e11 * .e58/(.e59^2 * (.e15)))) *
      .e1 * .e11), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * .e111 - 0.5 * .e30) * .e35/(.e15) + S * .e61 *
    .e40 * (.e35 * .e1/.e2 + .e35) * .e1 * epsilon1/((.e15) * .e20)^2) * epsilon1 -
    0.5 * (.e51/(.e15)))/(.e15) - .e61 * .e40 * .e35 * .e1/((.e15) * .e20)^2)/.e59 -
    2 * (.e51 * .e11 * .e63/(.e59^2 * (.e15)))) * .e2 * .e11), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e51 * .e64 * (.e16) * .e11 * .e12 * sqrt(.e16)/(.e65^2 *
    (.e15) * sqrt(.e15))))), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e51 * .e3 * .e11 * .e12 *
    .e69 * sqrt(.e16)/(.e70^2 * (.e15) * sqrt(.e15))))), uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e51 * .e4 * .e11 *
    .e12 * .e74 * sqrt(.e16)/(.e70^2 * (.e15) * sqrt(.e15))))), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e51 * .e75 *
    (.e17) * .e11 * .e13 * sqrt(.e17)/(.e76^2 * (.e15) * sqrt(.e15))))), Xvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e51 * .e5 *
    .e11 * .e13 * .e80 * sqrt(.e17)/(.e81^2 * (.e15) * sqrt(.e15))))), uHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e51 * .e6 *
    .e11 * .e13 * .e85 * sqrt(.e17)/(.e81^2 * (.e15) * sqrt(.e15))))), vHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e51 * .e86 *
    (.e18) * .e11 * .e14 * sqrt(.e18)/(.e87^2 * (.e15) * sqrt(.e15))))), Xvar)
  hessll[1:nXvar, (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e51 * .e7 *
    .e11 * .e14 * .e91 * sqrt(.e18)/(.e92^2 * (.e15) * sqrt(.e15))))), uHvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e51 * .e8 *
    .e11 * .e14 * .e96 * sqrt(.e18)/(.e92^2 * (.e15) * sqrt(.e15))))), vHvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e49 * .e51 *
    .e97 * (.e19) * .e11 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 + .e14) *
    (.e15) * sqrt(.e15))))), Xvar)
  hessll[1:nXvar, (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e49 * .e51 *
    .e9 * .e11 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 +
    .e13 + .e14) * (.e15) * sqrt(.e15))))), uHvar)
  hessll[1:nXvar, (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e49 * .e51 *
    .e10 * .e11 * .e105 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 +
    .e13 + .e14) * (.e15) * sqrt(.e15))))), vHvar)
  hessll[1:nXvar, (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * S * ((2 - 2 *
    (.e11/(1 + .e11 + .e12 + .e13 + .e14)))/.e106 - 2 * (.e107 * .e11/.e106^2)) *
    .e51 * .e11/((.e15) * sqrt(.e15)), Zvar)
  hessll[1:nXvar, (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e112 * .e51 * .e11 * .e12/((.e15) * sqrt(.e15)))), Zvar)
  hessll[1:nXvar, (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar *
    (-(S * .e113 * .e51 * .e11 * .e13/((.e15) * sqrt(.e15)))), Zvar)
  hessll[1:nXvar, (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Xvar * wHvar *
    (-(S * .e114 * .e51 * .e11 * .e14/((.e15) * sqrt(.e15)))), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e1 * (S * (0.5 * (S * .e35 * (S * (0.5 * (S * .e30 * epsilon1/(.e15)^2) -
    .e56 * .e40) * epsilon1 - 2 * (.e30/(.e15))) * epsilon1/(.e15)^2) - (0.5 *
    (S^2 * .e56 * .e35 * epsilon1^2/(.e15)^2) - (((0.5 * (.e1/(.e15)) + 1 - 0.5 *
    (0.5 * (1 - .e1/(.e15)) + .e1/(.e15))) * (1 - .e1/(.e15)) * .e2/.e20 + (2 -
    2 * (.e55^2 * .e1 * (.e15)/((.e15) * .e20)^2)) * .e20)/((.e15) * .e20)^2 +
    S^2 * .e56^2 * .e1 * epsilon1^2/((.e15) * .e20)) * .e35) * .e40) * epsilon1 -
    0.5 * ((S * .e57 * epsilon1 - .e35 * .e30/(.e15))/(.e15))) + S * .e57 * epsilon1 -
    0.5 * (.e35 * .e30/(.e15)))/.e59 - (0.5 * (.e53 * (1 + .e11 + .e12 + .e13 +
    .e14)/sqrt(.e15)) + 2 * (.e11 * .e58)) * .e1 * .e58/.e59^2) * .e1 * .e11),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 * ((1 - .e1/(.e15)) *
    .e2) - S^2 * .e61 * .e56 * .e1 * epsilon1^2)/(.e15) + 0.5 * ((.e1/(.e15) -
    1) * .e2/(.e15) + 1 - 0.5 * ((1 - .e1/(.e15)) * (1 - .e2/(.e15))))) * .e35/.e20 +
    0.5 * (S^2 * .e61 * .e35 * epsilon1^2/(.e15)^2)) * .e1 + .e61 * (1 - 2 *
    (.e55 * .e1 * (.e15) * .e20/((.e15) * .e20)^2)) * .e35) * .e40/((.e15) *
    .e20)^2 + 0.5 * (S * .e35 * (S * (0.5 * (S * .e30 * epsilon1/(.e15)^2) -
    .e56 * .e40) * epsilon1 - 2 * (.e30/(.e15))) * epsilon1/(.e15)^2)) * epsilon1 -
    0.5 * ((S * .e57 * epsilon1 - .e35 * .e30/(.e15))/(.e15)))/.e59 - (0.5 *
    (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e15)) + 2 * (.e11 * .e58)) *
    .e63/.e59^2) * .e1 * .e2 * .e11), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e64 * .e1 *
    (.e16) * .e11 * .e12 * .e58 * sqrt(.e16)/(.e65^2 * sqrt(.e15))))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e1 *
    .e3 * .e11 * .e12 * .e58 * .e69 * sqrt(.e16)/(.e70^2 * sqrt(.e15))))), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e1 *
    .e4 * .e11 * .e12 * .e74 * .e58 * sqrt(.e16)/(.e70^2 * sqrt(.e15))))), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (S * .e75 * .e1 * (.e17) * .e11 * .e13 * .e58 * sqrt(.e17)/(.e76^2 *
      sqrt(.e15))))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e1 * .e5 * .e11 * .e13 * .e58 * .e80 * sqrt(.e17)/(.e81^2 * sqrt(.e15))))),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar +
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e1 * .e6 * .e11 * .e13 * .e85 * .e58 * sqrt(.e17)/(.e81^2 * sqrt(.e15))))),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (S * .e86 * .e1 * (.e18) * .e11 * .e14 * .e58 * sqrt(.e18)/(.e87^2 *
      sqrt(.e15))))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e1 * .e7 * .e11 * .e14 * .e58 * .e91 * sqrt(.e18)/(.e92^2 * sqrt(.e15))))),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e1 * .e8 * .e11 * .e14 * .e96 * .e58 * sqrt(.e18)/(.e92^2 * sqrt(.e15))))),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    1):(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (S * .e49 * .e97 * .e1 * (.e19) * .e11 * .e58 * sqrt(.e19)/(.e98^2 *
      (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar +
    1):(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e49 * .e1 * .e9 * .e11 * .e58 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 *
      (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))))), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar)] <- crossprod(uHvar * wHvar *
    (-(4 * (.e49 * .e1 * .e10 * .e11 * .e105 * .e58 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 *
      (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))))), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(uHvar *
    wHvar * ((2 - 2 * (.e11/(1 + .e11 + .e12 + .e13 + .e14)))/.e106 - 2 * (.e107 *
    .e11/.e106^2)) * .e1 * .e11 * .e58/sqrt(.e15), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar *
    wHvar * (-(.e112 * .e1 * .e11 * .e12 * .e58/sqrt(.e15))), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    2 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar *
    wHvar * (-(.e113 * .e1 * .e11 * .e13 * .e58/sqrt(.e15))), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    3 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(uHvar *
    wHvar * (-(.e114 * .e1 * .e11 * .e14 * .e58/sqrt(.e15))), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e2/(.e15)) - 0.5 * (0.5 * (1 - .e2/(.e15)) + .e2/(.e15))) * (1 - .e2/(.e15)) +
    S^2 * .e61^2 * .e1 * .e2 * epsilon1^2/(((.e15) * .e20)^2 * (.e15))) * .e35 *
    .e1/.e20 + ((0.5 * (S^2 * .e35 * epsilon1^2/(.e15)^2) - 2 * (.e61 * .e35 *
    (.e15) * .e20/((.e15) * .e20)^2)) * .e2 + .e35) * .e61) * .e40 * .e1/((.e15) *
    .e20)^2 + S * (0.5 * (.e2 * (S * (.e61 * .e40 * .e1/((.e15) * .e20)^2 + 0.5 *
    (S * .e30 * epsilon1/(.e15)^2)) * epsilon1 - 2 * (.e30/(.e15)))) + 0.5 *
    .e30) * .e35 * epsilon1/(.e15)^2) * epsilon1 - (0.5 * (.e35 * .e30) + 0.5 *
    (.e2 * (S * .e62 * epsilon1 - .e35 * .e30/(.e15))))/(.e15))/.e59 - (0.5 *
    (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e15)) + 2 * (.e11 * .e63)) *
    .e2 * .e63/.e59^2) * .e2 * .e11), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (S * .e64 * (.e16) * .e2 * .e11 * .e12 * .e63 * sqrt(.e16)/(.e65^2 *
      sqrt(.e15))))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (.e3 * .e2 * .e11 * .e12 * .e63 * .e69 * sqrt(.e16)/(.e70^2 * sqrt(.e15))))),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e2 * .e4 * .e11 * .e12 * .e63 * .e74 * sqrt(.e16)/(.e70^2 *
    sqrt(.e15))))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (S * .e75 * (.e17) * .e2 * .e11 * .e13 * .e63 * sqrt(.e17)/(.e76^2 *
    sqrt(.e15))))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e5 * .e2 * .e11 * .e13 * .e63 * .e80 * sqrt(.e17)/(.e81^2 *
    sqrt(.e15))))), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e2 * .e6 * .e11 * .e13 * .e63 * .e85 * sqrt(.e17)/(.e81^2 *
    sqrt(.e15))))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (S * .e86 * (.e18) * .e2 * .e11 * .e14 * .e63 * sqrt(.e18)/(.e87^2 *
    sqrt(.e15))))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e7 * .e2 * .e11 * .e14 * .e63 * .e91 * sqrt(.e18)/(.e92^2 *
    sqrt(.e15))))), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e2 * .e8 * .e11 * .e14 * .e63 * .e96 * sqrt(.e18)/(.e92^2 *
    sqrt(.e15))))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (S * .e49 * .e97 * (.e19) * .e2 * .e11 * .e63 * sqrt(.e19)/(.e98^2 *
    (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))))), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e49 * .e9 * .e2 * .e11 * .e63 * .e102 * sqrt(.e19)/((.e53 *
    sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))))), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e49 * .e2 * .e10 * .e11 * .e63 * .e105 * sqrt(.e19)/((.e53 *
    sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 + .e14) * sqrt(.e15))))), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(vHvar *
    wHvar * ((2 - 2 * (.e11/(1 + .e11 + .e12 + .e13 + .e14)))/.e106 - 2 * (.e107 *
    .e11/.e106^2)) * .e2 * .e11 * .e63/sqrt(.e15), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar *
    wHvar * (-(.e112 * .e2 * .e11 * .e12 * .e63/sqrt(.e15))), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + 2 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 *
    nZHvar)] <- crossprod(vHvar * wHvar * (-(.e113 * .e2 * .e11 * .e13 * .e63/sqrt(.e15))),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + 3 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 *
    nZHvar)] <- crossprod(vHvar * wHvar * (-(.e114 * .e2 * .e11 * .e14 * .e63/sqrt(.e15))),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (((.e36 * (S * (.e41 * .e3/.e21 + S * .e31 * epsilon2) * epsilon2/(.e16) -
    .e31) + S * .e41 * (.e36 * .e3/.e4 + .e36) * .e3 * epsilon2/((.e16) * .e21))/.e65 -
    2 * (.e64^2 * .e12/.e65^2)) * .e12), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * ((.e67 * .e41 * .e36 + (S * ((0.5 * (S * (.e41 * .e3/.e21 +
    S * .e31 * epsilon2) * epsilon2/(.e16) - .e31) - 0.5 * .e31) * .e36/(.e16) -
    S * .e67 * .e41 * (.e36 * .e3/.e4 + .e36) * epsilon2) * epsilon2 - 0.5 *
    (.e64/(.e16)))/(.e16))/.e70 - 2 * (.e64 * .e12 * .e69/(.e70^2 * (.e16)))) *
    .e3 * .e12), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * (S * (.e41 * .e3/.e21 + S * .e31 * epsilon2) *
    epsilon2/(.e16) - .e31) - 0.5 * .e31) * .e36/(.e16) + S * .e72 * .e41 * (.e36 *
    .e3/.e4 + .e36) * .e3 * epsilon2/((.e16) * .e21)^2) * epsilon2 - 0.5 * (.e64/(.e16)))/(.e16) -
    .e72 * .e41 * .e36 * .e3/((.e16) * .e21)^2)/.e70 - 2 * (.e64 * .e12 * .e74/(.e70^2 *
    (.e16)))) * .e4 * .e12), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e64 * .e75 * (.e17) * .e12 * .e13 * sqrt(.e17)/(.e76^2 *
    (.e16) * sqrt(.e16))))), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e64 * .e5 * .e12 * .e13 * .e80 * sqrt(.e17)/(.e81^2 *
    (.e16) * sqrt(.e16))))), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e64 * .e6 * .e12 * .e13 * .e85 * sqrt(.e17)/(.e81^2 *
    (.e16) * sqrt(.e16))))), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (3 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e64 * .e86 * (.e18) * .e12 * .e14 * sqrt(.e18)/(.e87^2 *
    (.e16) * sqrt(.e16))))), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e64 * .e7 * .e12 * .e14 * .e91 * sqrt(.e18)/(.e92^2 *
    (.e16) * sqrt(.e16))))), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e64 * .e8 * .e12 * .e14 * .e96 * sqrt(.e18)/(.e92^2 *
    (.e16) * sqrt(.e16))))), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e49 * .e64 * .e97 * (.e19) * .e12 * sqrt(.e19)/(.e98^2 *
    (1 + .e11 + .e12 + .e13 + .e14) * (.e16) * sqrt(.e16))))), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (5 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e49 * .e64 * .e9 * .e12 * .e102 * sqrt(.e19)/((.e53 *
    sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 + .e14) * (.e16) * sqrt(.e16))))),
    uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (5 *
    nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (S * .e49 * .e64 * .e10 * .e12 * .e105 * sqrt(.e19)/((.e53 *
    sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 + .e14) * (.e16) * sqrt(.e16))))),
    vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    nZHvar)] <- crossprod(Xvar * wHvar * (-(S * (.e115 * .e64 * .e11 * .e12/((.e16) *
    sqrt(.e16))))), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * S * .e116 * .e64 *
    .e12/((.e16) * sqrt(.e16)), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S * .e113 * .e64 *
    .e12 * .e13/((.e16) * sqrt(.e16)))), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + 4 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S * .e114 * .e64 *
    .e12 * .e14/((.e16) * sqrt(.e16)))), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e3 * (S * (0.5 * (S * .e36 * (S * (0.5 * (S * .e31 * epsilon2/(.e16)^2) -
    .e67 * .e41) * epsilon2 - 2 * (.e31/(.e16))) * epsilon2/(.e16)^2) - (0.5 *
    (S^2 * .e67 * .e36 * epsilon2^2/(.e16)^2) - (((0.5 * (.e3/(.e16)) + 1 - 0.5 *
    (0.5 * (1 - .e3/(.e16)) + .e3/(.e16))) * (1 - .e3/(.e16)) * .e4/.e21 + (2 -
    2 * (.e66^2 * .e3 * (.e16)/((.e16) * .e21)^2)) * .e21)/((.e16) * .e21)^2 +
    S^2 * .e67^2 * .e3 * epsilon2^2/((.e16) * .e21)) * .e36) * .e41) * epsilon2 -
    0.5 * ((S * .e68 * epsilon2 - .e36 * .e31/(.e16))/(.e16))) + S * .e68 * epsilon2 -
    0.5 * (.e36 * .e31/(.e16)))/.e70 - (0.5 * (.e53 * (1 + .e11 + .e12 + .e13 +
    .e14)/sqrt(.e16)) + 2 * (.e12 * .e69)) * .e3 * .e69/.e70^2) * .e3 * .e12),
    uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((S * (((((0.5 * ((1 - .e3/(.e16)) * .e4) - S^2 * .e72 * .e67 *
    .e3 * epsilon2^2)/(.e16) + 0.5 * ((.e3/(.e16) - 1) * .e4/(.e16) + 1 - 0.5 *
    ((1 - .e3/(.e16)) * (1 - .e4/(.e16))))) * .e36/.e21 + 0.5 * (S^2 * .e72 *
    .e36 * epsilon2^2/(.e16)^2)) * .e3 + .e72 * (1 - 2 * (.e66 * .e3 * (.e16) *
    .e21/((.e16) * .e21)^2)) * .e36) * .e41/((.e16) * .e21)^2 + 0.5 * (S * .e36 *
    (S * (0.5 * (S * .e31 * epsilon2/(.e16)^2) - .e67 * .e41) * epsilon2 - 2 *
      (.e31/(.e16))) * epsilon2/(.e16)^2)) * epsilon2 - 0.5 * ((S * .e68 *
    epsilon2 - .e36 * .e31/(.e16))/(.e16)))/.e70 - (0.5 * (.e53 * (1 + .e11 +
    .e12 + .e13 + .e14)/sqrt(.e16)) + 2 * (.e12 * .e69)) * .e74/.e70^2) * .e3 *
    .e4 * .e12), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e75 * .e3 * (.e17) *
    .e12 * .e13 * .e69 * sqrt(.e17)/(.e76^2 * sqrt(.e16))))), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e3 * .e5 * .e12 * .e13 *
    .e69 * .e80 * sqrt(.e17)/(.e81^2 * sqrt(.e16))))), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e3 * .e6 * .e12 * .e13 *
    .e85 * .e69 * sqrt(.e17)/(.e81^2 * sqrt(.e16))))), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e86 * .e3 * (.e18) *
    .e12 * .e14 * .e69 * sqrt(.e18)/(.e87^2 * sqrt(.e16))))), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e3 * .e7 * .e12 * .e14 *
    .e69 * .e91 * sqrt(.e18)/(.e92^2 * sqrt(.e16))))), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e3 * .e8 * .e12 * .e14 *
    .e96 * .e69 * sqrt(.e18)/(.e92^2 * sqrt(.e16))))), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar + 4 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e49 * .e97 * .e3 *
    (.e19) * .e12 * .e69 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 + .e14) *
    sqrt(.e16))))), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 4 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e49 * .e3 * .e9 * .e12 *
    .e69 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e16))))), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 5 *
      nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e49 * .e3 * .e10 * .e12 *
    .e105 * .e69 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e16))))), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar + 5 *
      nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-((.e115 * .e3 * .e11 *
    .e12 * .e69/sqrt(.e16)))), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar + 5 * nuZUvar +
      5 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * .e116 * .e3 *
    .e12 * .e69/sqrt(.e16), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 * nXvar + 5 *
      nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e113 *
    .e3 * .e12 * .e13 * .e69/sqrt(.e16))), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 * nXvar + 5 *
      nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e114 *
    .e3 * .e12 * .e14 * .e69/sqrt(.e16))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 * (.e4/(.e16)) -
    0.5 * (0.5 * (1 - .e4/(.e16)) + .e4/(.e16))) * (1 - .e4/(.e16)) + S^2 * .e72^2 *
    .e3 * .e4 * epsilon2^2/(((.e16) * .e21)^2 * (.e16))) * .e36 * .e3/.e21 +
    ((0.5 * (S^2 * .e36 * epsilon2^2/(.e16)^2) - 2 * (.e72 * .e36 * (.e16) *
      .e21/((.e16) * .e21)^2)) * .e4 + .e36) * .e72) * .e41 * .e3/((.e16) *
    .e21)^2 + S * (0.5 * (.e4 * (S * (.e72 * .e41 * .e3/((.e16) * .e21)^2 + 0.5 *
    (S * .e31 * epsilon2/(.e16)^2)) * epsilon2 - 2 * (.e31/(.e16)))) + 0.5 *
    .e31) * .e36 * epsilon2/(.e16)^2) * epsilon2 - (0.5 * (.e36 * .e31) + 0.5 *
    (.e4 * (S * .e73 * epsilon2 - .e36 * .e31/(.e16))))/(.e16))/.e70 - (0.5 *
    (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e16)) + 2 * (.e12 * .e74)) *
    .e4 * .e74/.e70^2) * .e4 * .e12), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e75 * (.e17) * .e4 *
    .e12 * .e13 * .e74 * sqrt(.e17)/(.e76^2 * sqrt(.e16))))), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e5 * .e4 * .e12 * .e13 *
    .e74 * .e80 * sqrt(.e17)/(.e81^2 * sqrt(.e16))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e4 * .e6 * .e12 * .e13 *
    .e74 * .e85 * sqrt(.e17)/(.e81^2 * sqrt(.e16))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e86 * (.e18) * .e4 *
    .e12 * .e14 * .e74 * sqrt(.e18)/(.e87^2 * sqrt(.e16))))), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e7 * .e4 * .e12 * .e14 *
    .e74 * .e91 * sqrt(.e18)/(.e92^2 * sqrt(.e16))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e4 * .e8 * .e12 * .e14 *
    .e74 * .e96 * sqrt(.e18)/(.e92^2 * sqrt(.e16))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e49 * .e97 * (.e19) *
    .e4 * .e12 * .e74 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 + .e14) *
    sqrt(.e16))))), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e49 * .e9 * .e4 * .e12 *
    .e74 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e16))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e49 * .e4 * .e10 * .e12 *
    .e74 * .e105 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e16))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-((.e115 * .e4 * .e11 *
    .e12 * .e74/sqrt(.e16)))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * .e116 *
    .e4 * .e12 * .e74/sqrt(.e16), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e113 *
    .e4 * .e12 * .e13 * .e74/sqrt(.e16))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e114 *
    .e4 * .e12 * .e14 * .e74/sqrt(.e16))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (((.e37 * (S * (.e42 *
    .e5/.e22 + S * .e32 * epsilon3) * epsilon3/(.e17) - .e32) + S * .e42 * (.e37 *
    .e5/.e6 + .e37) * .e5 * epsilon3/((.e17) * .e22))/.e76 - 2 * (.e75^2 * .e13/.e76^2)) *
    .e13), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * ((.e78 * .e42 *
    .e37 + (S * ((0.5 * (S * (.e42 * .e5/.e22 + S * .e32 * epsilon3) * epsilon3/(.e17) -
    .e32) - 0.5 * .e32) * .e37/(.e17) - S * .e78 * .e42 * (.e37 * .e5/.e6 + .e37) *
    epsilon3) * epsilon3 - 0.5 * (.e75/(.e17)))/(.e17))/.e81 - 2 * (.e75 * .e13 *
    .e80/(.e81^2 * (.e17)))) * .e5 * .e13), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * (((S * ((0.5 *
    (S * (.e42 * .e5/.e22 + S * .e32 * epsilon3) * epsilon3/(.e17) - .e32) -
    0.5 * .e32) * .e37/(.e17) + S * .e83 * .e42 * (.e37 * .e5/.e6 + .e37) * .e5 *
    epsilon3/((.e17) * .e22)^2) * epsilon3 - 0.5 * (.e75/(.e17)))/(.e17) - .e83 *
    .e42 * .e37 * .e5/((.e17) * .e22)^2)/.e81 - 2 * (.e75 * .e13 * .e85/(.e81^2 *
    (.e17)))) * .e6 * .e13), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e75 * .e86 *
    (.e18) * .e13 * .e14 * sqrt(.e18)/(.e87^2 * (.e17) * sqrt(.e17))))), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e75 * .e7 *
    .e13 * .e14 * .e91 * sqrt(.e18)/(.e92^2 * (.e17) * sqrt(.e17))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e75 * .e8 *
    .e13 * .e14 * .e96 * sqrt(.e18)/(.e92^2 * (.e17) * sqrt(.e17))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e49 * .e75 *
    .e97 * (.e19) * .e13 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 + .e14) *
    (.e17) * sqrt(.e17))))), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e49 * .e75 *
    .e9 * .e13 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 +
    .e13 + .e14) * (.e17) * sqrt(.e17))))), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e49 * .e75 *
    .e10 * .e13 * .e105 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 +
    .e13 + .e14) * (.e17) * sqrt(.e17))))), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * (-(S * (.e115 *
    .e75 * .e11 * .e13/((.e17) * sqrt(.e17))))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e112 * .e75 * .e12 * .e13/((.e17) * sqrt(.e17)))), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar *
    S * .e117 * .e75 * .e13/((.e17) * sqrt(.e17)), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Xvar * wHvar *
    -(S * .e114 * .e75 * .e13 * .e14/((.e17) * sqrt(.e17))), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((.e5 * (S * (0.5 *
    (S * .e37 * (S * (0.5 * (S * .e32 * epsilon3/(.e17)^2) - .e78 * .e42) * epsilon3 -
      2 * (.e32/(.e17))) * epsilon3/(.e17)^2) - (0.5 * (S^2 * .e78 * .e37 *
    epsilon3^2/(.e17)^2) - (((0.5 * (.e5/(.e17)) + 1 - 0.5 * (0.5 * (1 - .e5/(.e17)) +
    .e5/(.e17))) * (1 - .e5/(.e17)) * .e6/.e22 + (2 - 2 * (.e77^2 * .e5 * (.e17)/((.e17) *
    .e22)^2)) * .e22)/((.e17) * .e22)^2 + S^2 * .e78^2 * .e5 * epsilon3^2/((.e17) *
    .e22)) * .e37) * .e42) * epsilon3 - 0.5 * ((S * .e79 * epsilon3 - .e37 *
    .e32/(.e17))/(.e17))) + S * .e79 * epsilon3 - 0.5 * (.e37 * .e32/(.e17)))/.e81 -
    (0.5 * (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e17)) + 2 * (.e13 *
      .e80)) * .e5 * .e80/.e81^2) * .e5 * .e13), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 *
    ((1 - .e5/(.e17)) * .e6) - S^2 * .e83 * .e78 * .e5 * epsilon3^2)/(.e17) +
    0.5 * ((.e5/(.e17) - 1) * .e6/(.e17) + 1 - 0.5 * ((1 - .e5/(.e17)) * (1 -
      .e6/(.e17))))) * .e37/.e22 + 0.5 * (S^2 * .e83 * .e37 * epsilon3^2/(.e17)^2)) *
    .e5 + .e83 * (1 - 2 * (.e77 * .e5 * (.e17) * .e22/((.e17) * .e22)^2)) * .e37) *
    .e42/((.e17) * .e22)^2 + 0.5 * (S * .e37 * (S * (0.5 * (S * .e32 * epsilon3/(.e17)^2) -
    .e78 * .e42) * epsilon3 - 2 * (.e32/(.e17))) * epsilon3/(.e17)^2)) * epsilon3 -
    0.5 * ((S * .e79 * epsilon3 - .e37 * .e32/(.e17))/(.e17)))/.e81 - (0.5 *
    (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e17)) + 2 * (.e13 * .e80)) *
    .e85/.e81^2) * .e5 * .e6 * .e13), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e86 * .e5 *
    (.e18) * .e13 * .e14 * .e80 * sqrt(.e18)/(.e87^2 * sqrt(.e17))))), Xvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e5 * .e7 *
    .e13 * .e14 * .e80 * .e91 * sqrt(.e18)/(.e92^2 * sqrt(.e17))))), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e5 * .e8 *
    .e13 * .e14 * .e96 * .e80 * sqrt(.e18)/(.e92^2 * sqrt(.e17))))), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e49 * .e97 *
    .e5 * (.e19) * .e13 * .e80 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e17))))), Xvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e49 * .e5 *
    .e9 * .e13 * .e80 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e17))))), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e49 * .e5 *
    .e10 * .e13 * .e105 * .e80 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e17))))), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-((.e115 *
    .e5 * .e11 * .e13 * .e80/sqrt(.e17)))), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e112 *
    .e5 * .e12 * .e13 * .e80/sqrt(.e17))), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar * wHvar *
    .e117 * .e5 * .e13 * .e80/sqrt(.e17), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(uHvar * wHvar *
    (-(.e114 * .e5 * .e13 * .e14 * .e80/sqrt(.e17))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e6/(.e17)) - 0.5 * (0.5 * (1 - .e6/(.e17)) + .e6/(.e17))) * (1 - .e6/(.e17)) +
    S^2 * .e83^2 * .e5 * .e6 * epsilon3^2/(((.e17) * .e22)^2 * (.e17))) * .e37 *
    .e5/.e22 + ((0.5 * (S^2 * .e37 * epsilon3^2/(.e17)^2) - 2 * (.e83 * .e37 *
    (.e17) * .e22/((.e17) * .e22)^2)) * .e6 + .e37) * .e83) * .e42 * .e5/((.e17) *
    .e22)^2 + S * (0.5 * (.e6 * (S * (.e83 * .e42 * .e5/((.e17) * .e22)^2 + 0.5 *
    (S * .e32 * epsilon3/(.e17)^2)) * epsilon3 - 2 * (.e32/(.e17)))) + 0.5 *
    .e32) * .e37 * epsilon3/(.e17)^2) * epsilon3 - (0.5 * (.e37 * .e32) + 0.5 *
    (.e6 * (S * .e84 * epsilon3 - .e37 * .e32/(.e17))))/(.e17))/.e81 - (0.5 *
    (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e17)) + 2 * (.e13 * .e85)) *
    .e6 * .e85/.e81^2) * .e6 * .e13), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e86 * (.e18) *
    .e6 * .e13 * .e14 * .e85 * sqrt(.e18)/(.e87^2 * sqrt(.e17))))), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e7 * .e6 *
    .e13 * .e14 * .e85 * .e91 * sqrt(.e18)/(.e92^2 * sqrt(.e17))))), uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e6 * .e8 *
    .e13 * .e14 * .e85 * .e96 * sqrt(.e18)/(.e92^2 * sqrt(.e17))))), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e49 * .e97 *
    (.e19) * .e6 * .e13 * .e85 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e17))))), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e49 * .e9 *
    .e6 * .e13 * .e85 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e17))))), uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e49 * .e6 *
    .e10 * .e13 * .e85 * .e105 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e17))))), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-((.e115 *
    .e6 * .e11 * .e13 * .e85/sqrt(.e17)))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e112 *
    .e6 * .e12 * .e13 * .e85/sqrt(.e17))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(vHvar * wHvar *
    .e117 * .e6 * .e13 * .e85/sqrt(.e17), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(vHvar * wHvar *
    (-(.e114 * .e6 * .e13 * .e14 * .e85/sqrt(.e17))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (((.e38 * (S * (.e43 *
    .e7/.e23 + S * .e33 * epsilon4) * epsilon4/(.e18) - .e33) + S * .e43 * (.e38 *
    .e7/.e8 + .e38) * .e7 * epsilon4/((.e18) * .e23))/.e87 - 2 * (.e86^2 * .e14/.e87^2)) *
    .e14), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * ((.e89 * .e43 *
    .e38 + (S * ((0.5 * (S * (.e43 * .e7/.e23 + S * .e33 * epsilon4) * epsilon4/(.e18) -
    .e33) - 0.5 * .e33) * .e38/(.e18) - S * .e89 * .e43 * (.e38 * .e7/.e8 + .e38) *
    epsilon4) * epsilon4 - 0.5 * (.e86/(.e18)))/(.e18))/.e92 - 2 * (.e86 * .e14 *
    .e91/(.e92^2 * (.e18)))) * .e7 * .e14), uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * (((S * ((0.5 *
    (S * (.e43 * .e7/.e23 + S * .e33 * epsilon4) * epsilon4/(.e18) - .e33) -
    0.5 * .e33) * .e38/(.e18) + S * .e94 * .e43 * (.e38 * .e7/.e8 + .e38) * .e7 *
    epsilon4/((.e18) * .e23)^2) * epsilon4 - 0.5 * (.e86/(.e18)))/(.e18) - .e94 *
    .e43 * .e38 * .e7/((.e18) * .e23)^2)/.e92 - 2 * (.e86 * .e14 * .e96/(.e92^2 *
    (.e18)))) * .e8 * .e14), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (.e49 * .e86 *
    .e97 * (.e19) * .e14 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 + .e14) *
    (.e18) * sqrt(.e18))))), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e49 * .e86 *
    .e9 * .e14 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 +
    .e13 + .e14) * (.e18) * sqrt(.e18))))), uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e49 * .e86 *
    .e10 * .e14 * .e105 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 + .e12 +
    .e13 + .e14) * (.e18) * sqrt(.e18))))), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * (-(S * (.e115 *
    .e86 * .e11 * .e14/((.e18) * sqrt(.e18))))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e112 * .e86 * .e12 * .e14/((.e18) * sqrt(.e18)))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar *
    (-(S * .e113 * .e86 * .e13 * .e14/((.e18) * sqrt(.e18)))), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Xvar * wHvar *
    (S * .e118 * .e86 * .e14/((.e18) * sqrt(.e18))), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 3 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((.e7 * (S * (0.5 *
    (S * .e38 * (S * (0.5 * (S * .e33 * epsilon4/(.e18)^2) - .e89 * .e43) * epsilon4 -
      2 * (.e33/(.e18))) * epsilon4/(.e18)^2) - (0.5 * (S^2 * .e89 * .e38 *
    epsilon4^2/(.e18)^2) - (((0.5 * (.e7/(.e18)) + 1 - 0.5 * (0.5 * (1 - .e7/(.e18)) +
    .e7/(.e18))) * (1 - .e7/(.e18)) * .e8/.e23 + (2 - 2 * (.e88^2 * .e7 * (.e18)/((.e18) *
    .e23)^2)) * .e23)/((.e18) * .e23)^2 + S^2 * .e89^2 * .e7 * epsilon4^2/((.e18) *
    .e23)) * .e38) * .e43) * epsilon4 - 0.5 * ((S * .e90 * epsilon4 - .e38 *
    .e33/(.e18))/(.e18))) + S * .e90 * epsilon4 - 0.5 * (.e38 * .e33/(.e18)))/.e92 -
    (0.5 * (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e18)) + 2 * (.e14 *
      .e91)) * .e7 * .e91/.e92^2) * .e7 * .e14), uHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 *
    ((1 - .e7/(.e18)) * .e8) - S^2 * .e94 * .e89 * .e7 * epsilon4^2)/(.e18) +
    0.5 * ((.e7/(.e18) - 1) * .e8/(.e18) + 1 - 0.5 * ((1 - .e7/(.e18)) * (1 -
      .e8/(.e18))))) * .e38/.e23 + 0.5 * (S^2 * .e94 * .e38 * epsilon4^2/(.e18)^2)) *
    .e7 + .e94 * (1 - 2 * (.e88 * .e7 * (.e18) * .e23/((.e18) * .e23)^2)) * .e38) *
    .e43/((.e18) * .e23)^2 + 0.5 * (S * .e38 * (S * (0.5 * (S * .e33 * epsilon4/(.e18)^2) -
    .e89 * .e43) * epsilon4 - 2 * (.e33/(.e18))) * epsilon4/(.e18)^2)) * epsilon4 -
    0.5 * ((S * .e90 * epsilon4 - .e38 * .e33/(.e18))/(.e18)))/.e92 - (0.5 *
    (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e18)) + 2 * (.e14 * .e91)) *
    .e96/.e92^2) * .e7 * .e8 * .e14), vHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e49 * .e97 *
    .e7 * (.e19) * .e14 * .e91 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e18))))), Xvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e49 * .e7 *
    .e9 * .e14 * .e91 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e18))))), uHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e49 * .e7 *
    .e10 * .e14 * .e105 * .e91 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e18))))), vHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-((.e115 *
    .e7 * .e11 * .e14 * .e91/sqrt(.e18)))), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e112 *
    .e7 * .e12 * .e14 * .e91/sqrt(.e18))), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar * wHvar *
    (-(.e113 * .e7 * .e13 * .e14 * .e91/sqrt(.e18))), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(uHvar * wHvar *
    .e118 * .e7 * .e14 * .e91/sqrt(.e18), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e8/(.e18)) - 0.5 * (0.5 * (1 - .e8/(.e18)) + .e8/(.e18))) * (1 - .e8/(.e18)) +
    S^2 * .e94^2 * .e7 * .e8 * epsilon4^2/(((.e18) * .e23)^2 * (.e18))) * .e38 *
    .e7/.e23 + ((0.5 * (S^2 * .e38 * epsilon4^2/(.e18)^2) - 2 * (.e94 * .e38 *
    (.e18) * .e23/((.e18) * .e23)^2)) * .e8 + .e38) * .e94) * .e43 * .e7/((.e18) *
    .e23)^2 + S * (0.5 * (.e8 * (S * (.e94 * .e43 * .e7/((.e18) * .e23)^2 + 0.5 *
    (S * .e33 * epsilon4/(.e18)^2)) * epsilon4 - 2 * (.e33/(.e18)))) + 0.5 *
    .e33) * .e38 * epsilon4/(.e18)^2) * epsilon4 - (0.5 * (.e38 * .e33) + 0.5 *
    (.e8 * (S * .e95 * epsilon4 - .e38 * .e33/(.e18))))/(.e18))/.e92 - (0.5 *
    (.e53 * (1 + .e11 + .e12 + .e13 + .e14)/sqrt(.e18)) + 2 * (.e14 * .e96)) *
    .e8 * .e96/.e92^2) * .e8 * .e14), vHvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (S * .e49 * .e97 *
    (.e19) * .e8 * .e14 * .e96 * sqrt(.e19)/(.e98^2 * (1 + .e11 + .e12 + .e13 +
    .e14) * sqrt(.e18))))), Xvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e49 * .e9 *
    .e8 * .e14 * .e96 * .e102 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e18))))), uHvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(vHvar * wHvar * (-(4 * (.e49 * .e8 *
    .e10 * .e14 * .e96 * .e105 * sqrt(.e19)/((.e53 * sqrt(.e19))^2 * (1 + .e11 +
    .e12 + .e13 + .e14) * sqrt(.e18))))), vHvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-((.e115 *
    .e8 * .e11 * .e14 * .e96/sqrt(.e18)))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e112 *
    .e8 * .e12 * .e14 * .e96/sqrt(.e18))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(vHvar * wHvar *
    (-(.e113 * .e8 * .e13 * .e14 * .e96/sqrt(.e18))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(vHvar * wHvar *
    .e118 * .e8 * .e14 * .e96/sqrt(.e18), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (((.e39 * (S * (.e44 *
    .e9/.e24 + S * .e34 * epsilon5) * epsilon5/(.e19) - .e34) + S * .e44 * (.e39 *
    .e9/.e10 + .e39) * .e9 * epsilon5/((.e19) * .e24))/.e98 - 2 * (.e49 * .e97^2/.e98^2)) *
    .e49), Xvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * ((.e100 * .e44 *
    .e39 + (S * ((0.5 * (S * (.e44 * .e9/.e24 + S * .e34 * epsilon5) * epsilon5/(.e19) -
    .e34) - 0.5 * .e34) * .e39/(.e19) - S * .e100 * .e44 * (.e39 * .e9/.e10 +
    .e39) * epsilon5) * epsilon5 - 0.5 * (.e97/(.e19)))/(.e19))/(.e53 * sqrt(.e19)) -
    2 * (.e49 * .e97 * .e102/((.e53 * sqrt(.e19))^2 * (.e19)))) * .e49 * .e9),
    uHvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(Xvar * wHvar * 2 * (S * (((S * ((0.5 *
    (S * (.e44 * .e9/.e24 + S * .e34 * epsilon5) * epsilon5/(.e19) - .e34) -
    0.5 * .e34) * .e39/(.e19) + S * .e103 * .e44 * (.e39 * .e9/.e10 + .e39) *
    .e9 * epsilon5/((.e19) * .e24)^2) * epsilon5 - 0.5 * (.e97/(.e19)))/(.e19) -
    .e103 * .e44 * .e39 * .e9/((.e19) * .e24)^2)/(.e53 * sqrt(.e19)) - 2 * (.e49 *
    .e97 * .e105/((.e53 * sqrt(.e19))^2 * (.e19)))) * .e49 * .e10), vHvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * (-(S * .e49 *
    .e119 * .e97 * .e11/((.e19) * sqrt(.e19)))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(Xvar * wHvar * (-(S *
    .e49 * .e120 * .e97 * .e12/((.e19) * sqrt(.e19)))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(Xvar * wHvar *
    (-(S * .e49 * .e122 * .e97 * .e13/((.e19) * sqrt(.e19)))), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Xvar * wHvar *
    (-(S * .e49 * .e121 * .e97 * .e14/((.e19) * sqrt(.e19)))), Zvar)
  hessll[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 4 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((.e9 * (S * (0.5 *
    (S * .e39 * (S * (0.5 * (S * .e34 * epsilon5/(.e19)^2) - .e100 * .e44) *
      epsilon5 - 2 * (.e34/(.e19))) * epsilon5/(.e19)^2) - (0.5 * (S^2 * .e100 *
    .e39 * epsilon5^2/(.e19)^2) - (((0.5 * (.e9/(.e19)) + 1 - 0.5 * (0.5 * (1 -
    .e9/(.e19)) + .e9/(.e19))) * (1 - .e9/(.e19)) * .e10/.e24 + (2 - 2 * (.e99^2 *
    .e9 * (.e19)/((.e19) * .e24)^2)) * .e24)/((.e19) * .e24)^2 + S^2 * .e100^2 *
    .e9 * epsilon5^2/((.e19) * .e24)) * .e39) * .e44) * epsilon5 - 0.5 * ((S *
    (0.5 * .e101 - .e100 * .e44 * .e39) * epsilon5 - .e39 * .e34/(.e19))/(.e19))) +
    S * (0.5 * .e101 - .e100 * .e44 * .e39) * epsilon5 - 0.5 * (.e39 * .e34/(.e19)))/(.e53 *
    sqrt(.e19)) - (0.5 * (.e53/sqrt(.e19)) + 2 * (.e49 * .e102)) * .e9 * .e102/(.e53 *
    sqrt(.e19))^2) * .e49 * .e9), uHvar)
  hessll[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 *
    ((1 - .e9/(.e19)) * .e10) - S^2 * .e103 * .e100 * .e9 * epsilon5^2)/(.e19) +
    0.5 * ((.e9/(.e19) - 1) * .e10/(.e19) + 1 - 0.5 * ((1 - .e9/(.e19)) * (1 -
      .e10/(.e19))))) * .e39/.e24 + 0.5 * (S^2 * .e103 * .e39 * epsilon5^2/(.e19)^2)) *
    .e9 + .e103 * (1 - 2 * (.e99 * .e9 * (.e19) * .e24/((.e19) * .e24)^2)) *
    .e39) * .e44/((.e19) * .e24)^2 + 0.5 * (S * .e39 * (S * (0.5 * (S * .e34 *
    epsilon5/(.e19)^2) - .e100 * .e44) * epsilon5 - 2 * (.e34/(.e19))) * epsilon5/(.e19)^2)) *
    epsilon5 - 0.5 * ((S * (0.5 * .e101 - .e100 * .e44 * .e39) * epsilon5 - .e39 *
    .e34/(.e19))/(.e19)))/(.e53 * sqrt(.e19)) - (0.5 * (.e53/sqrt(.e19)) + 2 *
    (.e49 * .e102)) * .e105/(.e53 * sqrt(.e19))^2) * .e49 * .e9 * .e10), vHvar)
  hessll[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-(.e49 * .e119 *
    .e9 * .e11 * .e102/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(uHvar * wHvar * (-(.e49 *
    .e120 * .e9 * .e12 * .e102/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(uHvar * wHvar *
    (-(.e49 * .e122 * .e9 * .e13 * .e102/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(uHvar * wHvar *
    (-(.e49 * .e121 * .e9 * .e14 * .e102/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e10/(.e19)) - 0.5 * (0.5 * (1 - .e10/(.e19)) + .e10/(.e19))) * (1 - .e10/(.e19)) +
    S^2 * .e103^2 * .e9 * .e10 * epsilon5^2/(((.e19) * .e24)^2 * (.e19))) * .e39 *
    .e9/.e24 + ((0.5 * (S^2 * .e39 * epsilon5^2/(.e19)^2) - 2 * (.e103 * .e39 *
    (.e19) * .e24/((.e19) * .e24)^2)) * .e10 + .e39) * .e103) * .e44 * .e9/((.e19) *
    .e24)^2 + S * (0.5 * (.e10 * (S * (.e103 * .e44 * .e9/((.e19) * .e24)^2 +
    0.5 * (S * .e34 * epsilon5/(.e19)^2)) * epsilon5 - 2 * (.e34/(.e19)))) +
    0.5 * .e34) * .e39 * epsilon5/(.e19)^2) * epsilon5 - (0.5 * (.e39 * .e34) +
    0.5 * (.e10 * (S * .e104 * epsilon5 - .e39 * .e34/(.e19))))/(.e19))/(.e53 *
    sqrt(.e19)) - (0.5 * (.e53/sqrt(.e19)) + 2 * (.e49 * .e105)) * .e10 * .e105/(.e53 *
    sqrt(.e19))^2) * .e49 * .e10), vHvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-(.e49 * .e119 *
    .e10 * .e11 * .e105/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(vHvar * wHvar * (-(.e49 *
    .e120 * .e10 * .e12 * .e105/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(vHvar * wHvar *
    (-(.e49 * .e122 * .e10 * .e13 * .e105/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(vHvar * wHvar *
    (-(.e49 * .e121 * .e10 * .e14 * .e105/sqrt(.e19))), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + nZHvar)] <- crossprod(Zvar * wHvar * ((1 - .e11/(1 +
    .e11 + .e12 + .e13 + .e14))/.e106 - 2 * (.e35 * .e11 * .e30/(.e106^2 * sqrt(.e15)))) *
    .e107 * .e11, Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e107/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2) + 2 * (.e108 *
    .e35 * .e30/(.e106^2 * sqrt(.e15)))) * .e11 * .e12)), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e107/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2) + 2 * (.e109 *
    .e35 * .e30/(.e106^2 * sqrt(.e15)))) * .e11 * .e13)), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar +
    1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e107/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2) + 2 * (.e110 *
    .e35 * .e30/(.e106^2 * sqrt(.e15)))) * .e11 * .e14)), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + 2 * nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)] <- crossprod(Zvar *
    wHvar * ((1 - .e12/(1 + .e11 + .e12 + .e13 + .e14))/.e106 - 2 * (.e36 * .e12 *
    .e31/(.e106^2 * sqrt(.e16)))) * .e108 * .e12, Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + 2 * nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    2 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e108/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2) + 2 * (.e109 *
    .e36 * .e31/(.e106^2 * sqrt(.e16)))) * .e12 * .e13)), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar + 5 *
    nuZUvar + 5 * nvZVvar + 2 * nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    3 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e108/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2) + 2 * (.e110 *
    .e36 * .e31/(.e106^2 * sqrt(.e16)))) * .e12 * .e14)), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    2 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)] <- crossprod(Zvar *
    wHvar * ((1 - .e13/(1 + .e11 + .e12 + .e13 + .e14))/.e106 - 2 * (.e37 * .e13 *
    .e32/(.e106^2 * sqrt(.e17)))) * .e109 * .e13, Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    3 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Zvar *
    wHvar * (-((.e109/(.e53 * (1 + .e11 + .e12 + .e13 + .e14)^2) + 2 * (.e110 *
    .e37 * .e32/(.e106^2 * sqrt(.e17)))) * .e13 * .e14)), Zvar)
  hessll[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar), (5 * nXvar + 5 * nuZUvar + 5 * nvZVvar +
    3 * nZHvar + 1):(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)] <- crossprod(Zvar *
    wHvar * ((1 - .e14/(1 + .e11 + .e12 + .e13 + .e14))/.e106 - 2 * (.e38 * .e14 *
    .e33/(.e106^2 * sqrt(.e18)))) * .e110 * .e14, Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for lcm 5 classes halfnormal-normal distribution
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
LCM5ChnormAlgOpt <- function(start, olsParam, dataTable, S, wHvar,
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csLCMfhalfnorm5C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(cLCMhalfnormlike5C(startVal, nXvar = nXvar,
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
  cat("LCM 5 Classes Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cLCMhalfnormlike5C(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike5C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cLCMhalfnormlike5C,
    grad = cgradLCMhalfnormlike5C, hess = chessLCMhalfnormlike5C,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cLCMhalfnormlike5C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessLCMhalfnormlike5C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradLCMhalfnormlike5C(mleObj$par,
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
      mleObj$hessian <- chessLCMhalfnormlike5C(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessLCMhalfnormlike5C(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cLCMhalfnormlike5C(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradLCMhalfnormlike5C(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) initHalf = initHalf))
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for lcmcross 5 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @param level level for confidence interval
#' @noRd
cLCM5Chalfnormeff <- function(object, level) {
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
  beta5 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar)]
  delta5 <- object$mlParam[(5 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    4 * object$nvZVvar)]
  phi5 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar)]
  theta1 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 2 * object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 3 * object$nZHvar)]
  theta4 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 3 * object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 4 * object$nZHvar)]
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
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon5 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta5), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  mustar5 <- -exp(Wu5) * object$S * epsilon5/(exp(Wu5) + exp(Wv5))
  sigmastar5 <- sqrt(exp(Wu5) * exp(Wv5)/(exp(Wu5) + exp(Wv5)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Pi5 <- 2/sqrt(exp(Wu5) + exp(Wv5)) * dnorm(object$S * epsilon5/sqrt(exp(Wu5) +
    exp(Wv5))) * pnorm(mustar5/sigmastar5)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc4 <- exp(Wz4)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc5 <- 1 - Probc1 - Probc2 - Probc3 - Probc4
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c5 <- Pi5 * Probc5/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4,
    Pcond_c5), 1, which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c ==
    2, Pcond_c2, ifelse(Group_c == 3, Pcond_c3, ifelse(Group_c ==
    4, Pcond_c4, Pcond_c5))))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
  u_c4 <- mustar4 + sigmastar4 * dnorm(mustar4/sigmastar4)/pnorm(mustar4/sigmastar4)
  u_c5 <- mustar5 + sigmastar5 * dnorm(mustar5/sigmastar5)/pnorm(mustar5/sigmastar5)
  u_c <- ifelse(Group_c == 1, u_c1, ifelse(Group_c == 2, u_c2,
    ifelse(Group_c == 3, u_c3, ifelse(Group_c == 4, u_c4,
      u_c5))))
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  ineff_c3 <- ifelse(Group_c == 3, u_c3, NA)
  ineff_c4 <- ifelse(Group_c == 4, u_c4, NA)
  ineff_c5 <- ifelse(Group_c == 5, u_c5, NA)
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
    teBC_c5 <- exp(-mustar5 + 1/2 * sigmastar5^2) * pnorm(mustar5/sigmastar5 -
      sigmastar5)/pnorm(mustar5/sigmastar5)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, ifelse(Group_c ==
      2, teBC_c2, ifelse(Group_c == 3, teBC_c3, ifelse(Group_c ==
      4, teBC_c4, teBC_c5))))
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    effBC_c3 <- ifelse(Group_c == 3, teBC_c3, NA)
    effBC_c4 <- ifelse(Group_c == 4, teBC_c4, NA)
    effBC_c5 <- ifelse(Group_c == 5, teBC_c5, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) *
      pnorm(mustar1/sigmastar1 + sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) *
      pnorm(mustar2/sigmastar2 + sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c3 <- exp(mustar3 + 1/2 * sigmastar3^2) *
      pnorm(mustar3/sigmastar3 + sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_reciprocal_c4 <- exp(mustar4 + 1/2 * sigmastar4^2) *
      pnorm(mustar4/sigmastar4 + sigmastar4)/pnorm(mustar4/sigmastar4)
    teBC_reciprocal_c5 <- exp(mustar5 + 1/2 * sigmastar5^2) *
      pnorm(mustar5/sigmastar5 + sigmastar5)/pnorm(mustar5/sigmastar5)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      ifelse(Group_c == 2, teBC_reciprocal_c2, ifelse(Group_c ==
        3, teBC_reciprocal_c3, ifelse(Group_c == 4, teBC_reciprocal_c4,
        teBC_reciprocal_c5))))
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    ReffBC_c3 <- ifelse(Group_c == 3, teBC_reciprocal_c3,
      NA)
    ReffBC_c4 <- ifelse(Group_c == 4, teBC_reciprocal_c4,
      NA)
    ReffBC_c5 <- ifelse(Group_c == 5, teBC_reciprocal_c5,
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
      teBC_reciprocal_c4 = teBC_reciprocal_c4, PosteriorProb_c5 = Pcond_c5,
      PriorProb_c5 = Probc5, u_c5 = u_c5, teBC_c5 = teBC_c5,
      teBC_reciprocal_c5 = teBC_reciprocal_c5, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, ineff_c4 = ineff_c4,
      ineff_c5 = ineff_c5, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      effBC_c3 = effBC_c3, effBC_c4 = effBC_c4, effBC_c5 = effBC_c5,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2, ReffBC_c3 = ReffBC_c3,
      ReffBC_c4 = ReffBC_c4, ReffBC_c5 = ReffBC_c5)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
      u_c3 = u_c3, PosteriorProb_c4 = Pcond_c4, PriorProb_c4 = Probc4,
      u_c4 = u_c4, PosteriorProb_c5 = Pcond_c5, PriorProb_c5 = Probc5,
      u_c5 = u_c5, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2,
      ineff_c3 = ineff_c3, ineff_c4 = ineff_c4, ineff_c5 = ineff_c5)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal effects for for lcmcross 5 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @noRd
cmargLCM5Chalfnorm_Eu <- function(object) {
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
  beta5 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar)]
  delta5 <- object$mlParam[(5 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    4 * object$nvZVvar)]
  phi5 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar)]
  theta1 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 2 * object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 3 * object$nZHvar)]
  theta4 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 3 * object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 4 * object$nZHvar)]
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
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon5 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta5), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  mustar5 <- -exp(Wu5) * object$S * epsilon5/(exp(Wu5) + exp(Wv5))
  sigmastar5 <- sqrt(exp(Wu5) * exp(Wv5)/(exp(Wu5) + exp(Wv5)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Pi5 <- 2/sqrt(exp(Wu5) + exp(Wv5)) * dnorm(object$S * epsilon5/sqrt(exp(Wu5) +
    exp(Wv5))) * pnorm(mustar5/sigmastar5)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc4 <- exp(Wz4)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc5 <- 1 - Probc1 - Probc2 - Probc3 - Probc4
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c5 <- Pi5 * Probc5/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4,
    Pcond_c5), 1, which.max)
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
  margEff_c5 <- kronecker(matrix(delta5[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu5/2) * dnorm(0), ncol = 1))
  colnames(margEff_c5) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c5")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      ifelse(Group_c == 2, margEff_c2[, c], ifelse(Group_c ==
        3, margEff_c3[, c], ifelse(Group_c == 4, margEff_c4[,
        c], margEff_c5[, c]))))
  }
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2,
    margEff_c3, margEff_c4, margEff_c5)
  return(margEff)
}

cmargLCM5Chalfnorm_Vu <- function(object) {
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
  beta5 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar)]
  delta5 <- object$mlParam[(5 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    4 * object$nvZVvar)]
  phi5 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    4 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar)]
  theta1 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 1):(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 2 * object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 3 * object$nZHvar)]
  theta4 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar +
    5 * object$nvZVvar + 3 * object$nZHvar + 1):(5 * object$nXvar +
    5 * object$nuZUvar + 5 * object$nvZVvar + 4 * object$nZHvar)]
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
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon5 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta5), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  mustar5 <- -exp(Wu5) * object$S * epsilon5/(exp(Wu5) + exp(Wv5))
  sigmastar5 <- sqrt(exp(Wu5) * exp(Wv5)/(exp(Wu5) + exp(Wv5)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) +
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) +
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Pi5 <- 2/sqrt(exp(Wu5) + exp(Wv5)) * dnorm(object$S * epsilon5/sqrt(exp(Wu5) +
    exp(Wv5))) * pnorm(mustar5/sigmastar5)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc4 <- exp(Wz4)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) +
    exp(Wz4))
  Probc5 <- 1 - Probc1 - Probc2 - Probc3 - Probc4
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Pcond_c5 <- Pi5 * Probc5/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 *
    Probc3 + Pi4 * Probc4 + Pi5 * Probc5)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4,
    Pcond_c5), 1, which.max)
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
  margEff_c5 <- kronecker(matrix(delta5[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu5) * (1 - (dnorm(0)/pnorm(0))^2),
    ncol = 1))
  colnames(margEff_c5) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c5")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      ifelse(Group_c == 2, margEff_c2[, c], ifelse(Group_c ==
        3, margEff_c3[, c], ifelse(Group_c == 4, margEff_c4[,
        c], margEff_c5[, c]))))
  }
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2,
    margEff_c3, margEff_c4, margEff_c5)
  return(margEff)
}
