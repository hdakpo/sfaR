################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: truncated normal - normal                                       #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for truncated normal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @noRd
ctruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- crossprod(matrix(omega), t(muHvar))[1, ]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  mustar <- (mu * exp(Wv) - exp(Wu) * S * epsilon)/(exp(Wu) +
    exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  ll <- (-1/2 * log(exp(Wu) + exp(Wv)) + dnorm((mu + S * epsilon)/sqrt(exp(Wu) +
    exp(Wv)), log = TRUE) + log(pnorm(mustar/sigmastar)) -
    log(pnorm(mu/sqrt(exp(Wu)))))
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for truncated normal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
csttruncnorm <- function(olsObj, epsiRes, S, nmuZUvar, nuZUvar,
  uHvar, muHvar, nvZVvar, vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  if (S * m3 > 0) {
    ## Coelli (1995) suggests 0.05 for gamma
    varu <- (abs(S * m3 * sqrt(pi/2)/(1 - 4/pi)))^(2/3)
  } else {
    varu <- (S * m3 * sqrt(pi/2)/(1 - 4/pi))^(2/3)
  }
  if (m2 < (pi - 2)/pi * varu) {
    varv <- abs(m2 - (1 - 2/pi) * varu)
  } else {
    varv <- m2 - (1 - 2/pi) * varu
  }
  dep_u <- 1/2 * log(((epsiRes^2 - varv) * pi/(pi - 2))^2)
  dep_v <- 1/2 * log((epsiRes^2 - (1 - 2/pi) * varu)^2)
  reg_hetu <- if (nuZUvar == 1) {
    lm(log(varu) ~ 1)
  } else {
    lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetu$coefficients)))
    stop("At least one of the OLS coefficients of 'uhet' is NA: ",
      paste(colnames(uHvar)[is.na(reg_hetu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  reg_hetv <- if (nvZVvar == 1) {
    lm(log(varv) ~ 1)
  } else {
    lm(dep_v ~ ., data = as.data.frame(vHvar[, 2:nvZVvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetv$coefficients)))
    stop("at least one of the OLS coefficients of 'vhet' is NA: ",
      paste(colnames(vHvar)[is.na(reg_hetv$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  reg_hetmu <- if (nmuZUvar == 1) {
    lm(epsiRes ~ 1)
  } else {
    lm(epsiRes ~ ., data = as.data.frame(muHvar[, 2:nmuZUvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetmu$coefficients)))
    stop("at least one of the OLS coefficients of 'muhet' is NA: ",
      paste(colnames(muHvar)[is.na(reg_hetmu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  omega <- coefficients(reg_hetmu)
  names(omega) <- paste0("Zmu_", colnames(muHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu * 2/pi), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, omega, delta, phi))
}

# Gradient of the likelihood function ----------
#' gradient for truncated normal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @noRd
cgradtruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- crossprod(matrix(omega), t(muHvar))[1, ]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
   .e1 <- Wu
  .e2 <- exp(.e1)
  .e3 <- exp(Wv)
  .e4 <- .e2 + .e3
  .e5 <- mu
  .e7 <- epsilon
  .e10 <- sqrt(.e2 * .e3/.e4)
  .e11 <- S * .e7
  .e12 <- .e4 * .e10
  .e13 <- .e5 + .e11
  .e15 <- .e5 * .e3 - S * .e2 * .e7
  .e16 <- .e15/.e12
  .e17 <- .e13/sqrt(.e4)
  .e19 <- exp(.e1/2)
  .e20 <- dnorm(.e16, 0, 1)
  .e21 <- dnorm(.e17)
  .e22 <- dnorm(.e17, 0, 1)
  .e23 <- .e5/.e19
  .e24 <- pnorm(.e16)
  .e25 <- .e12^2
  .e26 <- (0.5 * (.e22 * .e13^2/(.e21 * .e4)) - 0.5)/.e4
  .e28 <- .e22 * .e13/.e21
  .e29 <- dnorm(.e23, 0, 1)
  .e30 <- .e19 * pnorm(.e23)
  .e31 <- .e24 * .e10
  gradll <- cbind(S * Xvar * ((.e20 * .e2/.e31 + .e28)/.e4), muHvar * ((.e20 *
    .e3/.e31 - .e28)/.e4 - .e29/.e30), uHvar * ((.e26 - ((0.5 * ((1 - .e2/.e4) *
    .e3/.e10) + .e10) * .e15/.e25 + .e11/.e12) * .e20/.e24) * .e2 + 0.5 * (.e5 *
    .e29/.e30)), vHvar * ((.e26 + .e20 * (.e5/.e12 - (0.5 * ((1 - .e3/.e4) *
    .e2/.e10) + .e10) * .e15/.e25)/.e24) * .e3))
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for truncated normal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @noRd
chesstruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- crossprod(matrix(omega), t(muHvar))[1, ]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- Wu
  .e2 <- exp(Wv)
  .e3 <- exp(.e1)
  .e4 <- .e3 + .e2
  .e5 <- mu
  .e7 <- epsilon
  .e10 <- sqrt(.e3 * .e2/.e4)
  .e11 <- .e4 * .e10
  .e12 <- S * .e7
  .e13 <- .e5 * .e2
  .e14 <- .e5 + .e12
  .e16 <- S * .e3 * .e7
  .e17 <- .e13 - .e16
  .e18 <- .e17/.e11
  .e19 <- .e14/sqrt(.e4)
  .e20 <- dnorm(.e19)
  .e21 <- pnorm(.e18)
  .e22 <- dnorm(.e19, 0, 1)
  .e23 <- dnorm(.e18, 0, 1)
  .e24 <- .e14^2
  .e25 <- .e11^2
  .e27 <- exp(.e1/2)
  .e28 <- .e2/.e4
  .e29 <- .e3/.e4
  .e30 <- 1 - .e28
  .e31 <- 1 - .e29
  .e32 <- .e30 * .e3
  .e33 <- .e31 * .e2
  .e36 <- 0.5 * (.e32/.e10) + .e10
  .e37 <- .e20 * .e4
  .e38 <- .e5/.e27
  .e41 <- 0.5 * (.e33/.e10) + .e10
  .e42 <- .e21 * .e10
  .e43 <- .e22 * .e24
  .e44 <- .e5/.e11
  .e47 <- .e44 - .e36 * .e17/.e25
  .e48 <- .e12/.e11
  .e51 <- .e41 * .e17/.e25 + .e48
  .e52 <- pnorm(.e38)
  .e53 <- .e37^2
  .e54 <- .e42^2
  .e55 <- .e22/.e20
  .e56 <- dnorm(.e38, 0, 1)
  .e57 <- .e27 * .e52
  .e58 <- 0.5 * (.e43/.e37)
  .e60 <- .e23 * .e3
  .e63 <- ((0.5 - 0.5 * .e55) * .e24/.e4 - 1) * .e22 * .e14/.e20
  .e67 <- .e57^2
  .e68 <- .e18 + .e23/.e21
  .e69 <- .e24/.e4
  .e71 <- 0.5 * ((0.5 * (.e24/(.e20 * .e4^3)) - (0.5 * (.e43/.e4) + .e20)/.e53) *
    .e22 * .e24) - (.e58 - 0.5)/.e4
  .e72 <- .e23 * .e47
  .e73 <- .e23 * .e2
  .e74 <- .e43/.e53
  .e76 <- .e51 * .e23
  .e78 <- .e4 * .e21 * .e10
  .e79 <- .e17 * .e47
  .e80 <- .e73/.e42
  .e82 <- .e27^3 * .e52
  .e83 <- .e5 * .e56
  .e85 <- mu^2
  .e89 <- .e63 - .e60/.e42
  .e90 <- .e63 + .e80
  .e92 <- .e51 * .e68 * .e47
  .e93 <- .e51 * .e17
  .e96 <- ((1 - .e55) * .e24/.e4 - 1) * .e22/.e20
  .e101 <- ((.e55 - 1) * .e24/.e4 + 1) * .e22/.e20
  .e103 <- (.e17/.e42 + .e60 * .e2/.e54) * .e23/.e4
  .e104 <- .e71/.e4
  .e105 <- (0.5 * (.e33 * .e21/.e11) - .e76 * .e10) * .e3
  .e115 <- .e36/.e25 + .e72/.e78
  .e120 <- 0.5 * (((.e69 - 2)/.e37 - .e74) * .e22 * .e14/.e4)
  .e121 <- 0.5 * (((2 - .e69)/.e37 + .e74) * .e22 * .e14/.e4)
  .e122 <- 0.5 * (.e31 * .e30)
  .e124 <- 0.5 * (.e32 * .e21/.e11) + .e72 * .e10
  .e126 <- 0.5 * (.e85/.e82) - (0.5 * .e57 - 0.5 * .e83)/.e67
  .e127 <- 1/.e10
  .e128 <- 2 * (.e41 * .e4 * .e17 * .e10/.e25)
  .e129 <- 2 * (.e36 * .e4 * .e17 * .e10/.e25)
  hessll <- matrix(0, nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar, ncol = nXvar +
    nmuZUvar + nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * ((.e96 - (.e17/(.e2 * .e21 * .e10) +
    .e60/.e54) * .e23 * .e3/.e4)/.e4) * wHvar, Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(S * Xvar * (-((.e96 +
    .e103)/.e4)) * wHvar, muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(S *
    Xvar * ((.e120 - ((.e41/.e25 - .e76/.e78) * .e3 - (.e93/.e2 + .e127)/.e4) *
    .e23/.e21) * .e3) * wHvar, uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(S * Xvar * ((.e120 - (.e115 * .e3 + .e79/(.e4 * .e2)) *
    .e23/.e21) * .e2) * wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(muHvar *
    (.e56 * (.e56/.e67 + .e5/.e82) - (.e101 + (.e17/(.e3 * .e21 * .e10) + .e73/.e54) *
      .e23 * .e2/.e4)/.e4) * wHvar, muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(muHvar * ((.e121 - (.e41 * .e2/.e25 - .e51 * (.e17/.e3 +
    .e80)/.e4) * .e23/.e21) * .e3 + 0.5 * (((1 - .e85/.e27^2)/.e57 - .e83/.e67) *
    .e56)) * wHvar, uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(muHvar * ((((.e127 - .e79/.e3)/.e4 -
    .e115 * .e2) * .e23/.e21 + .e121) * .e2) * wHvar, vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(uHvar * (((.e71 * .e3 + .e58 -
    0.5)/.e4 - ((.e41 * (.e13 - (.e128 + 3 * .e12) * .e3) + (0.5 * .e29 - 0.5 *
    (0.5 * .e31 + .e29)) * .e31 * .e2 * .e17/.e10)/.e25 + .e51^2 * .e68 * .e3 +
    .e48) * .e23/.e21) * .e3 + 0.5 * (.e5 * .e126 * .e56)) * wHvar, uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(uHvar *
    (((.e92 - ((0.5 * (.e33/.e4) + 0.5 * ((.e29 - 1) * .e2/.e4 + 1 - .e122)) *
      .e17/.e10 + .e5 * .e41 - .e36 * (.e128 + .e12))/.e25) * .e23/.e21 + .e104) *
      .e3 * .e2) * wHvar, vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(vHvar *
    (((.e71 * .e2 + .e58 - 0.5)/.e4 + .e23 * (.e44 - ((((3 * .e5 - .e129) * .e2 -
      .e16) * .e36 + (0.5 * .e28 - 0.5 * (0.5 * .e30 + .e28)) * .e30 * .e3 *
      .e17/.e10)/.e25 + .e68 * .e2 * .e47^2))/.e21) * .e2) * wHvar, vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for truncated normal-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param method algorithm for solver
#' @param derivs type of derivatives
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @param accuracy accuracy for numerical derivatives
#' @param stepsize stepsize for numerical derivatives
#' @noRd
truncnormAlgOpt <- function(start, randStart, sdStart, olsParam, dataTable, S, nXvar,
  muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar, method,
  derivs, printInfo, itermax, stepmax, tol, gradtol, hessianType, qac, accuracy,
  stepsize) {
   ## starting values and log likelihood ------
  startVal <- if (!is.null(start))
    start else csttruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], S = S,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar)
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- ctruncnormlike(startVal, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar)
  ## automatic differentiation ------
  if (derivs == "ad") {
    if (method %in% c("bfgs", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
    }
    cgradbhhhtruncnormAD <- RTMB::MakeADFun(function(p) ctruncnormlike(p, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar),
    startVal, ADreport = TRUE, silent = TRUE)
  cgradhesstruncnormAD <- RTMB::MakeADFun(function(p) ctruncnormlike(p, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar),
    startVal, ADreport = FALSE, silent = TRUE)
    fnGradObstrunc <- function(parm) {
      Ta1 <- RTMB::GetTape(cgradbhhhtruncnormAD)
      Ta2 <- RTMB::MakeTape(function(weight) {
        WT <- RTMB::MakeTape(function(x) sum(Ta1(x) * weight), parm)
        (WT$jacfun())(RTMB::advector(parm))
      }, rep(1, nrow(Xvar)))
      t((Ta2$jacfun())(RTMB::advector(rep(1, nrow(Xvar)))))
    }
    ### solve for different algorithms ------
    mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(p) -cgradhesstruncnormAD$fn(p),
      gr = function(p) -cgradhesstruncnormAD$gr(p), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
        maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
      maxLikAlgo = maxRoutine(fn = cgradhesstruncnormAD$fn, grad = cgradhesstruncnormAD$gr,
        hess = cgradhesstruncnormAD$he, start = startVal, finalHessian = TRUE,
        control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax,
          reltol = tol, tol = tol, qac = qac)), bhhh = maxLik::maxBHHH(fn = cgradbhhhtruncnormAD$fn,
        grad = fnGradObstrunc, hess = cgradhesstruncnormAD$he, start = startVal,
        finalHessian = TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac)), sr1 = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesstruncnormAD$fn(p), gr = function(p) -cgradhesstruncnormAD$gr(p),
        method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesstruncnormAD$fn(p), gr = function(p) -cgradhesstruncnormAD$gr(p),
        hs = function(p) as(-cgradhesstruncnormAD$he(p), "dgCMatrix"), method = "Sparse",
        control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(p) -cgradhesstruncnormAD$fn(p),
        gr = function(p) -cgradhesstruncnormAD$gr(p), hess = function(p) -cgradhesstruncnormAD$he(p),
        print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
      nlminb = nlminb(start = startVal, objective = function(p) -cgradhesstruncnormAD$fn(p),
        gradient = function(p) -cgradhesstruncnormAD$gr(p), hessian = function(p) -cgradhesstruncnormAD$he(p),
        control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
          eval.max = itermax, rel.tol = tol, x.tol = tol)))
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$gradient <- cgradhesstruncnormAD$gr(mleObj$par)
    }
    mlParam <- if (method %in% c("ucminf", "nlminb")) {
      mleObj$par
    } else {
      if (method %in% c("maxLikAlgo", "bhhh")) {
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
        mleObj$hessian <- cgradhesstruncnormAD$he(mleObj$par)
      if (method == "sr1")
        mleObj$hessian <- cgradhesstruncnormAD$he(mleObj$solution)
      if (method == "mla")
        mleObj$hessian <- cgradhesstruncnormAD$he(mleObj$b)
    }
    mleObj$logL_OBS <- cgradbhhhtruncnormAD$fn(mlParam)
    mleObj$gradL_OBS <- fnGradObstrunc(mlParam)
    rm(cgradbhhhtruncnormAD, cgradhesstruncnormAD, fnGradObstrunc)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
      maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
        bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
        nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
        sann = function(...) maxLik::maxSANN(...))
      method <- "maxLikAlgo"
    }
    ## numerical derivatives ------
    if (derivs == "numerical") {
      cgradtruncnormlikeNum <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
        muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar) {
        calculus::jacobian(ctruncnormlike, var = unname(parm), params = list(nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar), accuracy = accuracy, stepsize = stepsize)
      }
      chesstruncnormlikeNum <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
        muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar) {
        calculus::hessian(ctruncnormlike, var = unname(parm), params = list(nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar), accuracy = accuracy, stepsize = stepsize)
      }
      ### solve for different algorithms ------
      mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ctruncnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar)), gr = function(parm) -colSums(cgradtruncnormlikeNum(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
        maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
        maxLikAlgo = maxRoutine(fn = ctruncnormlike, grad = cgradtruncnormlikeNum,
          hess = chesstruncnormlikeNum, start = startVal, finalHessian = if (hessianType ==
          2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ctruncnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar)), gr = function(parm) -colSums(cgradtruncnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar)), method = "SR1", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)),
          gr = function(parm) -colSums(cgradtruncnormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar)), hs = function(parm) as(-chesstruncnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
          fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)),
          gr = function(parm) -colSums(cgradtruncnormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar)), hess = function(parm) -chesstruncnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar), print.info = printInfo, maxiter = itermax,
          epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar)), gradient = function(parm) -colSums(cgradtruncnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar)), hessian = function(parm) -chesstruncnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar), control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
          eval.max = itermax, rel.tol = tol, x.tol = tol)))
      if (method %in% c("ucminf", "nlminb")) {
        mleObj$gradient <- colSums(cgradtruncnormlikeNum(mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar))
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
          mleObj$hessian <- chesstruncnormlikeNum(parm = mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar)
        if (method == "sr1")
          mleObj$hessian <- chesstruncnormlikeNum(parm = mleObj$solution,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar)
      }
      mleObj$logL_OBS <- ctruncnormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
      mleObj$gradL_OBS <- cgradtruncnormlikeNum(parm = mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
    } else {
      ## analytical derivatives ------
      if (derivs == "analytical") {
        ### solve for different algorithms ------
        mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
          fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)),
          gr = function(parm) -colSums(cgradtruncnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
          maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
          maxLikAlgo = maxRoutine(fn = ctruncnormlike, grad = cgradtruncnormlike,
          hess = chesstruncnormlike, start = startVal, finalHessian = if (hessianType ==
            2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
            iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ctruncnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar)), gr = function(parm) -colSums(cgradtruncnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
          muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
          S = S, wHvar = wHvar)), method = "SR1", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
            nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)),
          gr = function(parm) -colSums(cgradtruncnormlike(parm, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar)), hs = function(parm) as(-chesstruncnormlike(parm,
            nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar), "dgCMatrix"), method = "Sparse",
          control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
            prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
            preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
          fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
            nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
            vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)),
          gr = function(parm) -colSums(cgradtruncnormlike(parm, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar)), hess = function(parm) -chesstruncnormlike(parm,
            nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar), print.info = printInfo,
          maxiter = itermax, epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar)), gradient = function(parm) -colSums(cgradtruncnormlike(parm,
            nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar)), hessian = function(parm) -chesstruncnormlike(parm,
            nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar), control = list(iter.max = itermax,
            trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
            x.tol = tol)))
        if (method %in% c("ucminf", "nlminb")) {
          mleObj$gradient <- colSums(cgradtruncnormlike(mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar))
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
          mleObj$hessian <- chesstruncnormlike(parm = mleObj$par, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar)
          if (method == "sr1")
          mleObj$hessian <- chesstruncnormlike(parm = mleObj$solution, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
            muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
            Xvar = Xvar, S = S, wHvar = wHvar)
        }
        mleObj$logL_OBS <- ctruncnormlike(mlParam, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar)
        mleObj$gradL_OBS <- cgradtruncnormlike(parm = mlParam, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
          wHvar = wHvar)
      }
    }
  }
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for truncated normal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
ctruncnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  mu <- crossprod(matrix(omega), t(muHvar))[1, ]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta), t(Xvar))[1, ]
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * epsilon)/(exp(Wu) +
    exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    m <- ifelse(mustar > 0, mustar, 0)
    teMO <- exp(-m)
    teBC <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBCLB <- exp(-uUB)
    teBCUB <- exp(-uLB)
    teBC_reciprocal <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    res <- data.frame(u = u, uLB = uLB, uUB = uUB, teJLMS = teJLMS,
      m = m, teMO = teMO, teBC = teBC, teBCLB = teBCLB,
      teBCUB = teBCUB, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u, uLB = uLB, uUB = uUB, m = m)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for truncated normal-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargtruncnorm_Eu <- function(object) {
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  mu <- crossprod(matrix(omega), t(muHvar))[1, ]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Lambda <- mu/exp(Wu/2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(1 - Lambda * dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2,
      ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2)/2 * ((1 + Lambda^2) * dnorm(Lambda)/pnorm(Lambda) +
      Lambda * (dnorm(Lambda)/pnorm(Lambda))^2), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  colnames(margEff) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(data.frame(margEff))
}

cmargtruncnorm_Vu <- function(object) {
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  mu <- crossprod(matrix(omega), t(muHvar))[1, ]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Lambda <- mu/exp(Wu/2)
  m1 <- exp(Wu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  m2 <- exp(Wu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) -
    (dnorm(Lambda)/pnorm(Lambda))^2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(1/exp(Wu/2) * dnorm(Lambda)/pnorm(Lambda) * (m1^2 -
      m2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - 1/2 * dnorm(Lambda)/pnorm(Lambda) *
      (Lambda + Lambda^3 + (2 + 3 * Lambda^2) * dnorm(Lambda)/pnorm(Lambda) +
        2 * Lambda * (dnorm(Lambda)/pnorm(Lambda))^2)),
      ncol = 1))
  idTRUE_mu <- (substring(names(omega)[-1], 5)) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  colnames(margEff) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(data.frame(margEff))
}
