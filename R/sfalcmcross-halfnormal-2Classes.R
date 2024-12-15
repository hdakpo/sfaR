################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Latent Class Stochastic Frontier Analysis                             #
# Number of Classes: 2L                                                        #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for lcm 2 classes halfnormal-normal distribution
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
cLCMhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- crossprod(matrix(delta1), t(uHvar))[1, ]
  Wu2 <- crossprod(matrix(delta2), t(uHvar))[1, ]
  Wv1 <- crossprod(matrix(phi1), t(vHvar))[1, ]
  Wv2 <- crossprod(matrix(phi2), t(vHvar))[1, ]
  Wz <- crossprod(matrix(theta), t(Zvar))[1, ]
  epsilon1 <- Yvar - crossprod(matrix(beta1), t(Xvar))[1, ]
  epsilon2 <- Yvar - crossprod(matrix(beta2), t(Xvar))[1, ]
  mustar1 <- -exp(Wu1) * S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ll <- log(Probc1 * Pi1 + Probc2 * Pi2)
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
  # ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for lcm 2 classes halfnormal-normal distribution
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
csLCMfhalfnorm2C <- function(olsObj, epsiRes, nXvar, nuZUvar,
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
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    paste0("Cl1_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for lcm 2 classes halfnormal-normal distribution
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
cgradLCMhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- crossprod(matrix(delta1), t(uHvar))[1, ]
  Wu2 <- crossprod(matrix(delta2), t(uHvar))[1, ]
  Wv1 <- crossprod(matrix(phi1), t(vHvar))[1, ]
  Wv2 <- crossprod(matrix(phi2), t(vHvar))[1, ]
  Wz <- crossprod(matrix(theta), t(Zvar))[1, ]
  epsilon1 <- Yvar - crossprod(matrix(beta1), t(Xvar))[1, ]
  epsilon2 <- Yvar - crossprod(matrix(beta2), t(Xvar))[1, ]
  .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wz)
  .e6 <- .e1 + .e2
  .e7 <- sqrt(.e6)
  .e8 <- sqrt(.e1 * .e2/(.e6))
  .e9 <- (S * .e1 * epsilon1/((.e6) * .e8))
  .e10 <- dnorm(S * epsilon1/.e7)
  .e11 <- dnorm(-.e9)
  .e12 <- pnorm(-.e9)
  .e13 <- .e3 + .e4
  .e14 <- sqrt(.e13)
  .e15 <- sqrt(.e3 * .e4/(.e13))
  .e16 <- dnorm(S * epsilon2/.e14)
  .e17 <- (S * .e3 * epsilon2/((.e13) * .e15))
  .e18 <- pnorm(-.e17)
  .e19 <- dnorm(-.e17)
  .e20 <- (.e11 * .e10 * .e1/.e8 + S * .e10 * .e12 * epsilon1)
  .e21 <- (1 - .e5/(1 + .e5))
  .e22 <- (.e21 * .e16 * .e18/.e14)
  .e23 <- (.e10 * .e5 * .e12/((1 + .e5) * .e7))
  .e24 <- ((1 + .e5) * (2 * .e22 + 2 * .e23) * (.e6) * .e7)
  .e25 <- ((2 * .e22 + 2 * .e23) * .e7)
  .e26 <- ((1 + .e5) * .e10 * .e12/((1 + .e5) * .e7)^2)
  .e27 <- ((1 - .e1/(.e6)) * .e2/.e8)
  .e28 <- (S * .e10 * .e12 * epsilon1/(.e6)^2)
  .e29 <- (1/((.e6) * .e8) - (0.5 * .e27 + .e8) * .e1/((.e6) * .e8)^2)
  .e30 <- (S * (0.5 * .e28 - .e29 * .e11 * .e10) * epsilon1/(1 + .e5) - 0.5 * .e26)
  .e31 <- ((1 - .e2/(.e6)) * .e1/.e8)
  .e32 <- ((0.5 * .e31 + .e8) * .e11 * .e10 * .e1/((.e6) * .e8)^2 + 0.5 * .e28)
  .e33 <- (S * .e32 * epsilon1/(1 + .e5) - 0.5 * .e26)
  .e34 <- (.e19 * .e16 * .e3/.e15 + S * .e16 * .e18 * epsilon2)
  .e35 <- ((2 * .e22 + 2 * .e23) * (.e13) * .e14)
  .e36 <- (S * .e16 * .e18 * epsilon2/(.e13)^2)
  .e37 <- (0.5 * ((1 - .e3/(.e13)) * .e4/.e15) + .e15)
  .e38 <- (1/((.e13) * .e15) - .e37 * .e3/((.e13) * .e15)^2)
  .e39 <- (0.5 * .e36 - .e38 * .e19 * .e16)
  .e40 <- (S * .e39 * epsilon2 - 0.5 * (.e16 * .e18/(.e13)))
  .e41 <- ((2 * .e22 + 2 * .e23) * .e14)
  .e42 <- ((1 - .e4/(.e13)) * .e3/.e15)
  .e43 <- ((0.5 * .e42 + .e15) * .e19 * .e16 * .e3/((.e13) * .e15)^2 + 0.5 * .e36)
  .e44 <- (S * .e43 * epsilon2 - 0.5 * (.e16 * .e18/(.e13)))
  .e45 <- (1/((1 + .e5) * .e7) - .e5 * .e7/((1 + .e5) * .e7)^2)
  .e46 <- (.e45 * .e10 * .e12)
  .e47 <- (.e21 * .e16 * .e18/((1 + .e5) * .e14))
  gradll <- cbind(2 * (S * Xvar * .e20 * .e5/.e24), 2 * (uHvar * .e1 * .e5 * .e30/.e25),
    2 * (vHvar * .e2 * .e5 * .e33/.e25), 2 * (S * Xvar * .e21 * .e34/.e35), 2 *
      (uHvar * .e21 * .e3 * .e40/.e41), 2 * (vHvar * .e21 * .e4 * .e44/.e41),
    Zvar * (2 * .e46 - 2 * .e47) * .e5/(2 * .e22 + 2 * .e23))
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for lcm 2 classes halfnormal-normal distribution
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
chessLCMhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- crossprod(matrix(delta1), t(uHvar))[1, ]
  Wu2 <- crossprod(matrix(delta2), t(uHvar))[1, ]
  Wv1 <- crossprod(matrix(phi1), t(vHvar))[1, ]
  Wv2 <- crossprod(matrix(phi2), t(vHvar))[1, ]
  Wz <- crossprod(matrix(theta), t(Zvar))[1, ]
  epsilon1 <- Yvar - crossprod(matrix(beta1), t(Xvar))[1, ]
  epsilon2 <- Yvar - crossprod(matrix(beta2), t(Xvar))[1, ]
   .e1 <- exp(Wu1)
  .e2 <- exp(Wv1)
  .e3 <- exp(Wu2)
  .e4 <- exp(Wv2)
  .e5 <- exp(Wz)
  .e6 <- .e1 + .e2
  .e7 <- sqrt(.e6)
  .e8 <- sqrt(.e1 * .e2/(.e6))
  .e9 <- (S * .e1 * epsilon1/((.e6) * .e8))
  .e10 <- dnorm(S * epsilon1/.e7)
  .e11 <- dnorm(-.e9)
  .e12 <- pnorm(-.e9)
  .e13 <- .e3 + .e4
  .e14 <- sqrt(.e13)
  .e15 <- sqrt(.e3 * .e4/(.e13))
  .e16 <- dnorm(S * epsilon2/.e14)
  .e17 <- (S * .e3 * epsilon2/((.e13) * .e15))
  .e18 <- pnorm(-.e17)
  .e19 <- dnorm(-.e17)
  .e20 <- (.e11 * .e10 * .e1/.e8 + S * .e10 * .e12 * epsilon1)
  .e21 <- (1 - .e5/(1 + .e5))
  .e22 <- (.e21 * .e16 * .e18/.e14)
  .e23 <- (.e10 * .e5 * .e12/((1 + .e5) * .e7))
  .e24 <- ((1 + .e5) * (2 * .e22 + 2 * .e23) * (.e6) * .e7)
  .e25 <- ((2 * .e22 + 2 * .e23) * .e7)
  .e26 <- ((1 + .e5) * .e10 * .e12/((1 + .e5) * .e7)^2)
  .e27 <- ((1 - .e1/(.e6)) * .e2/.e8)
  .e28 <- (S * .e10 * .e12 * epsilon1/(.e6)^2)
  .e29 <- (1/((.e6) * .e8) - (0.5 * .e27 + .e8) * .e1/((.e6) * .e8)^2)
  .e30 <- (S * (0.5 * .e28 - .e29 * .e11 * .e10) * epsilon1/(1 + .e5) - 0.5 * .e26)
  .e31 <- ((1 - .e2/(.e6)) * .e1/.e8)
  .e32 <- ((0.5 * .e31 + .e8) * .e11 * .e10 * .e1/((.e6) * .e8)^2 + 0.5 * .e28)
  .e33 <- (S * .e32 * epsilon1/(1 + .e5) - 0.5 * .e26)
  .e34 <- (.e19 * .e16 * .e3/.e15 + S * .e16 * .e18 * epsilon2)
  .e35 <- ((2 * .e22 + 2 * .e23) * (.e13) * .e14)
  .e36 <- (S * .e16 * .e18 * epsilon2/(.e13)^2)
  .e37 <- (0.5 * ((1 - .e3/(.e13)) * .e4/.e15) + .e15)
  .e38 <- (1/((.e13) * .e15) - .e37 * .e3/((.e13) * .e15)^2)
  .e39 <- (0.5 * .e36 - .e38 * .e19 * .e16)
  .e40 <- (S * .e39 * epsilon2 - 0.5 * (.e16 * .e18/(.e13)))
  .e41 <- ((2 * .e22 + 2 * .e23) * .e14)
  .e42 <- ((1 - .e4/(.e13)) * .e3/.e15)
  .e43 <- ((0.5 * .e42 + .e15) * .e19 * .e16 * .e3/((.e13) * .e15)^2 + 0.5 * .e36)
  .e44 <- (S * .e43 * epsilon2 - 0.5 * (.e16 * .e18/(.e13)))
  .e45 <- (1/((1 + .e5) * .e7) - .e5 * .e7/((1 + .e5) * .e7)^2)
  .e46 <- (.e45 * .e10 * .e12)
  .e47 <- (.e21 * .e16 * .e18/((1 + .e5) * .e14))
  .e48 <- (S * (.e11 * .e1/.e8 + S * .e12 * epsilon1) * epsilon1/(.e6) - .e12)
  .e49 <- (.e10 * .e1/.e2 + .e10)
  .e50 <- ((1 + .e5) * .e20/(((1 + .e5) * .e7)^2 * (.e6)))
  .e51 <- (.e41^2 * (1 + .e5) * (.e6) * .e7)
  .e52 <- (S * (0.5 * .e28 - .e29 * .e11 * .e10) * epsilon1 - (1 + .e5)^2 * .e10 *
    .e12/((1 + .e5) * .e7)^2)
  .e53 <- ((1 + .e5) * .e52/((1 + .e5) * .e7)^2)
  .e54 <- (0.5 * (S * .e12 * epsilon1/(.e6)^2) - .e29 * .e11)
  .e55 <- (S * .e54 * epsilon1 - 2 * (.e12/(.e6)))
  .e56 <- (S * .e10 * .e55 * epsilon1/(.e6)^2)
  hessll <- matrix(nrow = 2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar, ncol = 2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * wHvar * 2 * (((.e10 * .e48 + S *
    .e11 * .e49 * .e1 * epsilon1/((.e6) * .e8))/.e24 - 2 * (.e20^2 * .e5/.e24^2)) *
    .e5), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(Xvar * wHvar * 2 *
    (S * (((.e29 * .e11 * .e10 + S * ((0.5 * .e48 - 0.5 * .e12) * .e10/(.e6) -
      S * .e29 * .e11 * .e49 * epsilon1) * epsilon1/(.e6))/(1 + .e5) - 0.5 *
      .e50)/.e25 - 2 * (.e20 * .e5 * .e30/(.e25^2 * (1 + .e5) * (.e6)))) *
      .e1 * .e5), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * .e48 - 0.5 * .e12) * .e10/(.e6) + S * (0.5 *
    .e31 + .e8) * .e11 * .e49 * .e1 * epsilon1/((.e6) * .e8)^2) * epsilon1/(.e6) -
    (0.5 * .e31 + .e8) * .e11 * .e10 * .e1/((.e6) * .e8)^2)/(1 + .e5) - 0.5 *
    .e50)/.e25 - 2 * (.e20 * .e5 * .e33/(.e25^2 * (1 + .e5) * (.e6)))) * .e2 *
    .e5), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * (-(4 * (.e21 * .e20 * .e34 * (.e13) * .e5 * .e14/(.e35^2 * (1 + .e5) *
    (.e6) * .e7)))), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e21 * .e20 * .e3 * .e5 *
    .e40 * .e14/.e51))), uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(Xvar * wHvar * (-(4 * (S * .e21 * .e20 * .e4 *
    .e5 * .e44 * .e14/.e51))), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(Xvar * wHvar * (S * (2 * .e45 -
    2 * ((2 * .e46 - 2 * .e47) * .e5/((1 + .e5) * (2 * .e22 + 2 * .e23) * .e7))) *
    .e20 * .e5/((2 * .e22 + 2 * .e23) * (.e6))), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e1 * (S * (0.5 * .e56 - (0.5 * (.e29 * .e10 * epsilon1^2/(.e6)^2) -
    (((0.5 * (.e1/(.e6)) + 1 - 0.5 * (0.5 * (1 - .e1/(.e6)) + .e1/(.e6))) * (1 -
      .e1/(.e6)) * .e2/.e8 + (2 - 2 * ((0.5 * .e27 + .e8)^2 * .e1 * (.e6)/((.e6) *
      .e8)^2)) * .e8)/((.e6) * .e8)^2 + .e29^2 * .e1 * epsilon1^2/((.e6) *
      .e8)) * .e10) * .e11) * epsilon1/(1 + .e5) - 0.5 * .e53) + S * (0.5 *
    .e28 - .e29 * .e11 * .e10) * epsilon1/(1 + .e5) - 0.5 * .e26)/.e25 - (0.5 *
    ((2 * .e22 + 2 * .e23)/.e7) + 2 * (.e5 * .e30)) * .e1 * .e30/.e25^2) * .e1 *
    .e5), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * wHvar * 2 * (((S * (((((0.5 * ((1 - .e1/(.e6)) *
    .e2) - (0.5 * .e31 + .e8) * .e29 * .e1 * epsilon1^2)/(.e6) + 0.5 * ((.e1/(.e6) -
    1) * .e2/(.e6) + 1 - 0.5 * ((1 - .e1/(.e6)) * (1 - .e2/(.e6))))) * .e10/.e8 +
    0.5 * ((0.5 * .e31 + .e8) * .e10 * epsilon1^2/(.e6)^2)) * .e1 + (0.5 * .e31 +
    .e8) * (1 - 2 * ((0.5 * .e27 + .e8) * .e1 * (.e6) * .e8/((.e6) * .e8)^2)) *
    .e10) * .e11/((.e6) * .e8)^2 + 0.5 * .e56) * epsilon1/(1 + .e5) - 0.5 * .e53)/.e25 -
    (0.5 * ((2 * .e22 + 2 * .e23)/.e7) + 2 * (.e5 * .e30)) * .e33/.e25^2) * .e1 *
    .e2 * .e5), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (S * .e21 * .e34 *
    .e1 * (.e13) * .e5 * .e30 * .e14/(.e35^2 * .e7)))), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e21 *
    .e1 * .e3 * .e5 * .e30 * .e40 * .e14/(.e41^2 * .e7)))), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar * wHvar * (-(4 * (.e21 *
    .e1 * .e4 * .e5 * .e44 * .e30 * .e14/(.e41^2 * .e7)))), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(uHvar *
    wHvar * (2 * ((0.5 * (.e45 * .e10 * epsilon1^2/(.e6)^2) - ((0.5/.e7 - (1 +
    .e5)^2 * .e7/((1 + .e5) * .e7)^2) * .e5 + 0.5 * ((1 + .e5)/.e7)) * .e10/((1 +
    .e5) * .e7)^2) * .e12 - S * .e45 * .e29 * .e11 * .e10 * epsilon1) - 2 * ((2 *
    .e46 - 2 * .e47) * .e5 * .e30/.e25)) * .e1 * .e5/(2 * .e22 + 2 * .e23), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 *
    (.e2/(.e6)) - 0.5 * (0.5 * (1 - .e2/(.e6)) + .e2/(.e6))) * (1 - .e2/(.e6)) +
    (0.5 * .e31 + .e8)^2 * .e1 * .e2 * epsilon1^2/(((.e6) * .e8)^2 * (.e6))) *
    .e10 * .e1/.e8 + ((0.5 * (.e10 * epsilon1^2/(.e6)^2) - 2 * ((0.5 * .e31 +
    .e8) * .e10 * (.e6) * .e8/((.e6) * .e8)^2)) * .e2 + .e10) * (0.5 * .e31 +
    .e8)) * .e11 * .e1/((.e6) * .e8)^2 + S * (0.5 * (.e2 * (S * ((0.5 * .e31 +
    .e8) * .e11 * .e1/((.e6) * .e8)^2 + 0.5 * (S * .e12 * epsilon1/(.e6)^2)) *
    epsilon1 - 2 * (.e12/(.e6)))) + 0.5 * .e12) * .e10 * epsilon1/(.e6)^2) *
    epsilon1/(1 + .e5) - (0.5 * (.e10 * .e12) + 0.5 * (.e2 * (S * .e32 * epsilon1 -
    (1 + .e5)^2 * .e10 * .e12/((1 + .e5) * .e7)^2))) * (1 + .e5)/((1 + .e5) *
    .e7)^2)/.e25 - (0.5 * ((2 * .e22 + 2 * .e23)/.e7) + 2 * (.e5 * .e33)) * .e2 *
    .e33/.e25^2) * .e2 * .e5), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (S * .e21 * .e34 * (.e13) * .e2 * .e5 * .e33 * .e14/(.e35^2 * .e7)))),
    Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(vHvar * wHvar *
    (-(4 * (.e21 * .e3 * .e2 * .e5 * .e33 * .e40 * .e14/(.e41^2 * .e7)))), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(vHvar *
    wHvar * (-(4 * (.e21 * .e2 * .e4 * .e5 * .e33 * .e44 * .e14/(.e41^2 * .e7)))),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(vHvar *
    wHvar * (2 * ((0.5 * (.e45 * .e10 * epsilon1^2/(.e6)^2) - ((0.5/.e7 - (1 +
    .e5)^2 * .e7/((1 + .e5) * .e7)^2) * .e5 + 0.5 * ((1 + .e5)/.e7)) * .e10/((1 +
    .e5) * .e7)^2) * .e12 + S * (0.5 * .e31 + .e8) * .e45 * .e11 * .e10 * .e1 *
    epsilon1/((.e6) * .e8)^2) - 2 * ((2 * .e46 - 2 * .e47) * .e5 * .e33/.e25)) *
    .e2 * .e5/(2 * .e22 + 2 * .e23), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (((.e16 * (S * (.e19 * .e3/.e15 + S * .e18 * epsilon2) * epsilon2/(.e13) -
    .e18) + S * .e19 * (.e16 * .e3/.e4 + .e16) * .e3 * epsilon2/((.e13) * .e15))/.e35 -
    2 * (.e21 * .e34^2/.e35^2)) * .e21), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * ((.e38 * .e19 * .e16 + (S * ((0.5 * (S * (.e19 * .e3/.e15 +
    S * .e18 * epsilon2) * epsilon2/(.e13) - .e18) - 0.5 * .e18) * .e16/(.e13) -
    S * .e38 * .e19 * (.e16 * .e3/.e4 + .e16) * epsilon2) * epsilon2 - 0.5 *
    (.e34/(.e13)))/(.e13))/.e41 - 2 * (.e21 * .e34 * .e40/(.e41^2 * (.e13)))) *
    .e21 * .e3), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(Xvar *
    wHvar * 2 * (S * (((S * ((0.5 * (S * (.e19 * .e3/.e15 + S * .e18 * epsilon2) *
    epsilon2/(.e13) - .e18) - 0.5 * .e18) * .e16/(.e13) + S * (0.5 * .e42 + .e15) *
    .e19 * (.e16 * .e3/.e4 + .e16) * .e3 * epsilon2/((.e13) * .e15)^2) * epsilon2 -
    0.5 * (.e34/(.e13)))/(.e13) - (0.5 * .e42 + .e15) * .e19 * .e16 * .e3/((.e13) *
    .e15)^2)/.e41 - 2 * (.e21 * .e34 * .e44/(.e41^2 * (.e13)))) * .e21 * .e4),
    vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar), (2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(Xvar * wHvar * (-(S * .e21 * (2 * ((2 * .e46 - 2 *
    .e47)/(2 * .e22 + 2 * .e23)) + 2/(1 + .e5)) * .e34 * .e5/.e35)), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((.e3 * (S * (0.5 * (S * .e16 * (S * (0.5 * (S * .e18 * epsilon2/(.e13)^2) -
    .e38 * .e19) * epsilon2 - 2 * (.e18/(.e13))) * epsilon2/(.e13)^2) - (0.5 *
    (.e38 * .e16 * epsilon2^2/(.e13)^2) - (((0.5 * (.e3/(.e13)) + 1 - 0.5 * (0.5 *
    (1 - .e3/(.e13)) + .e3/(.e13))) * (1 - .e3/(.e13)) * .e4/.e15 + (2 - 2 *
    (.e37^2 * .e3 * (.e13)/((.e13) * .e15)^2)) * .e15)/((.e13) * .e15)^2 + .e38^2 *
    .e3 * epsilon2^2/((.e13) * .e15)) * .e16) * .e19) * epsilon2 - 0.5 * ((S *
    .e39 * epsilon2 - .e16 * .e18/(.e13))/(.e13))) + S * .e39 * epsilon2 - 0.5 *
    (.e16 * .e18/(.e13)))/.e41 - (0.5 * ((2 * .e22 + 2 * .e23)/.e14) + 2 * (.e21 *
    .e40)) * .e3 * .e40/.e41^2) * .e21 * .e3), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(uHvar *
    wHvar * 2 * (((S * (((((0.5 * ((1 - .e3/(.e13)) * .e4) - (0.5 * .e42 + .e15) *
    .e38 * .e3 * epsilon2^2)/(.e13) + 0.5 * ((.e3/(.e13) - 1) * .e4/(.e13) +
    1 - 0.5 * ((1 - .e3/(.e13)) * (1 - .e4/(.e13))))) * .e16/.e15 + 0.5 * ((0.5 *
    .e42 + .e15) * .e16 * epsilon2^2/(.e13)^2)) * .e3 + (0.5 * .e42 + .e15) *
    (1 - 2 * (.e37 * .e3 * (.e13) * .e15/((.e13) * .e15)^2)) * .e16) * .e19/((.e13) *
    .e15)^2 + 0.5 * (S * .e16 * (S * (0.5 * (S * .e18 * epsilon2/(.e13)^2) -
    .e38 * .e19) * epsilon2 - 2 * (.e18/(.e13))) * epsilon2/(.e13)^2)) * epsilon2 -
    0.5 * ((S * .e39 * epsilon2 - .e16 * .e18/(.e13))/(.e13)))/.e41 - (0.5 *
    ((2 * .e22 + 2 * .e23)/.e14) + 2 * (.e21 * .e40)) * .e44/.e41^2) * .e21 *
    .e3 * .e4), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar),
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
      nvZVvar + nZHvar)] <- crossprod(uHvar * wHvar * (-(.e21 * (2 * ((2 *
    .e46 - 2 * .e47) * .e40/(2 * .e22 + 2 * .e23)) + 2 * (S * .e39 * epsilon2/(1 +
    .e5) - 0.5 * ((1 + .e5) * .e16 * .e18/((1 + .e5) * .e14)^2))) * .e3 * .e5/.e41)),
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(vHvar * wHvar * 2 * (((S * ((((0.5 * (.e4/(.e13)) -
    0.5 * (0.5 * (1 - .e4/(.e13)) + .e4/(.e13))) * (1 - .e4/(.e13)) + (0.5 *
    .e42 + .e15)^2 * .e3 * .e4 * epsilon2^2/(((.e13) * .e15)^2 * (.e13))) * .e16 *
    .e3/.e15 + ((0.5 * (.e16 * epsilon2^2/(.e13)^2) - 2 * ((0.5 * .e42 + .e15) *
    .e16 * (.e13) * .e15/((.e13) * .e15)^2)) * .e4 + .e16) * (0.5 * .e42 + .e15)) *
    .e19 * .e3/((.e13) * .e15)^2 + S * (0.5 * (.e4 * (S * ((0.5 * .e42 + .e15) *
    .e19 * .e3/((.e13) * .e15)^2 + 0.5 * (S * .e18 * epsilon2/(.e13)^2)) * epsilon2 -
    2 * (.e18/(.e13)))) + 0.5 * .e18) * .e16 * epsilon2/(.e13)^2) * epsilon2 -
    (0.5 * (.e16 * .e18) + 0.5 * (.e4 * (S * .e43 * epsilon2 - .e16 * .e18/(.e13))))/(.e13))/.e41 -
    (0.5 * ((2 * .e22 + 2 * .e23)/.e14) + 2 * (.e21 * .e44)) * .e4 * .e44/.e41^2) *
    .e21 * .e4), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 *
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(vHvar * wHvar * (-(.e21 * (2 * ((2 *
    .e46 - 2 * .e47) * .e44/(2 * .e22 + 2 * .e23)) + 2 * (S * .e43 * epsilon2/(1 +
    .e5) - 0.5 * ((1 + .e5) * .e16 * .e18/((1 + .e5) * .e14)^2))) * .e4 * .e5/.e41)),
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(Zvar * wHvar * ((2 * (.e21 *
    (1/((1 + .e5)^2 * .e14) + .e14/((1 + .e5) * .e14)^2) * .e16 * .e18) - ((2 *
    .e46 - 2 * .e47)^2/(2 * .e22 + 2 * .e23) + 2 * ((2 - 2 * ((1 + .e5) * (.e6) *
    .e5/((1 + .e5) * .e7)^2)) * .e10 * .e12 * .e7/((1 + .e5) * .e7)^2))) * .e5 +
    2 * .e46 - 2 * .e47) * .e5/(2 * .e22 + 2 * .e23), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for lcm 2 classes halfnormal-normal distribution
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
LCM2ChnormAlgOpt <- function(start, olsParam, dataTable, S, wHvar,
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar,
  Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csLCMfhalfnorm2C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(cLCMhalfnormlike2C(startVal, nXvar = nXvar,
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
  cat("LCM 2 Classes Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cLCMhalfnormlike2C(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cLCMhalfnormlike2C,
    grad = cgradLCMhalfnormlike2C, hess = chessLCMhalfnormlike2C,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cLCMhalfnormlike2C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessLCMhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradLCMhalfnormlike2C(mleObj$par,
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
      mleObj$hessian <- chessLCMhalfnormlike2C(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessLCMhalfnormlike2C(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cLCMhalfnormlike2C(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradLCMhalfnormlike2C(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) initHalf = initHalf))
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for lcmcross 2 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @param level level for confidence interval
#' @noRd
cLCM2Chalfnormeff <- function(object, level) {
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
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- crossprod(matrix(delta1), t(uHvar))[1, ]
  Wu2 <- crossprod(matrix(delta2), t(uHvar))[1, ]
  Wv1 <- crossprod(matrix(phi1), t(vHvar))[1, ]
  Wv2 <- crossprod(matrix(phi2), t(vHvar))[1, ]
  Wz <- crossprod(matrix(theta), t(Zvar))[1, ]
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta1), t(Xvar))[1, ]
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta2), t(Xvar))[1, ]
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
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
    teJLMS_c <- exp(-u_c)
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
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal effects for for lcmcross 2 classes halfnormal-normal distribution
#' @param object object of class lcmcross
#' @noRd
cmargLCM2Chalfnorm_Eu <- function(object) {
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
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- crossprod(matrix(delta1), t(uHvar))[1, ]
  Wu2 <- crossprod(matrix(delta2), t(uHvar))[1, ]
  Wv1 <- crossprod(matrix(phi1), t(vHvar))[1, ]
  Wv2 <- crossprod(matrix(phi2), t(vHvar))[1, ]
  Wz <- crossprod(matrix(theta), t(Zvar))[1, ]
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta1), t(Xvar))[1, ]
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta2), t(Xvar))[1, ]
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu1/2) * dnorm(0), ncol = 1))
  colnames(margEff_c1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu2/2) * dnorm(0), ncol = 1))
  colnames(margEff_c2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      margEff_c2[, c])
  }
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2)
  return(margEff)
}

cmargLCM2Chalfnorm_Vu <- function(object) {
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
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- crossprod(matrix(delta1), t(uHvar))[1, ]
  Wu2 <- crossprod(matrix(delta2), t(uHvar))[1, ]
  Wv1 <- crossprod(matrix(phi1), t(vHvar))[1, ]
  Wv2 <- crossprod(matrix(phi2), t(vHvar))[1, ]
  Wz <- crossprod(matrix(theta), t(Zvar))[1, ]
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta1), t(Xvar))[1, ]
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta2), t(Xvar))[1, ]
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
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
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in seq_len(ncol(margEff_c1))) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c],
      margEff_c2[, c])
  }
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2)
  return(margEff)
}
