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
  .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- .e1 + .e2
  .e4 <- sqrt(.e3)
  .e5 <- sqrt(.e1 * .e2/(.e3))
  .e6 <- ((.e3) * .e5)
  .e7 <- (mu * .e2 - S * .e1 * epsilon)
  .e8 <- dnorm(.e7/.e6)
  .e9 <- pnorm(.e7/.e6)
  .e10 <- (mu + S * epsilon)
  .e11 <- exp(Wu/2)
  .e12 <- dnorm(mu/.e11)
  .e13 <- pnorm(mu/.e11)
  .e14 <- dnorm(.e10/.e4)
  .e15 <- (0.5 * ((1 - .e2/(.e3)) * .e1/.e5) + .e5)
  .e16 <- (mu/.e6 - .e15 * .e7/.e6^2)
  .e17 <- (.e14 * .e10^2/(.e14 * (.e3)))
  .e18 <- ((1 - .e1/(.e3)) * .e2/.e5)
  .e19 <- ((0.5 * .e18 + .e5) * .e7/.e6^2 + S * epsilon/.e6)
  gradll <- cbind(S * Xvar * (.e8 * .e1/(.e9 * .e5) + .e14 * .e10/.e14)/(.e3),
    muHvar * ((.e8 * .e2/(.e9 * .e5) - .e14 * .e10/.e14)/(.e3) - .e12/(.e11 *
      .e13)), uHvar * (((0.5 * .e17 - 0.5)/(.e3) - .e19 * .e8/.e9) * .e1 +
      0.5 * (mu * .e12/(.e11 * .e13))), vHvar * ((0.5 * .e17 - 0.5)/(.e3) +
      .e8 * .e16/.e9) * .e2)
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
  .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- .e1 + .e2
  .e4 <- sqrt(.e3)
  .e5 <- sqrt(.e1 * .e2/(.e3))
  .e6 <- ((.e3) * .e5)
  .e7 <- (mu * .e2 - S * .e1 * epsilon)
  .e8 <- dnorm(.e7/.e6)
  .e9 <- pnorm(.e7/.e6)
  .e10 <- (mu + S * epsilon)
  .e11 <- exp(Wu/2)
  .e12 <- dnorm(mu/.e11)
  .e13 <- pnorm(mu/.e11)
  .e14 <- dnorm(.e10/.e4)
  .e15 <- (0.5 * ((1 - .e2/(.e3)) * .e1/.e5) + .e5)
  .e16 <- (mu/.e6 - .e15 * .e7/.e6^2)
  .e17 <- (.e14 * .e10^2/(.e14 * (.e3)))
  .e18 <- ((1 - .e1/(.e3)) * .e2/.e5)
  .e19 <- ((0.5 * .e18 + .e5) * .e7/.e6^2 + S * epsilon/.e6)
  .e20 <- ((1 - .e14/.e14) * .e10^2/(.e3) - 1)
  .e21 <- ((.e10^2/(.e3) - 2)/(.e14 * (.e3)) - .e14 * .e10^2/(.e14 * (.e3))^2)
  .e22 <- (.e21 * .e14 * .e10/(.e3))
  .e23 <- (.e15/.e6^2 + .e8 * .e16/((.e3) * .e9 * .e5))
  .e24 <- (((2 - .e10^2/(.e3))/(.e14 * (.e3)) + .e14 * .e10^2/(.e14 * (.e3))^2) *
    .e14 * .e10/(.e3))
  .e25 <- (0.5 * (.e14 * .e10^2/(.e3)) + .e14)
  .e26 <- (.e10^2/(.e14 * (.e3)^3))
  .e27 <- (0.5 * .e26 - .e25/(.e14 * (.e3))^2)
  .e28 <- (0.5 * (.e27 * .e14 * .e10^2) - (0.5 * .e17 - 0.5)/(.e3))
  .e29 <- ((0.5 * .e18 + .e5) * (.e3) * .e7 * .e5/.e6^2)
  hessll <- matrix(0, nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar, ncol = nXvar +
    nmuZUvar + nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * ((.e20 * .e14/.e14 - (.e7/(.e2 *
    .e9 * .e5) + .e8 * .e1/(.e9 * .e5)^2) * .e8 * .e1/(.e3))/(.e3)) * wHvar,
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(S * Xvar * (-((.e20 *
    .e14/.e14 + (.e7/(.e9 * .e5) + .e8 * .e1 * .e2/(.e9 * .e5)^2) * .e8/(.e3))/(.e3))) *
    wHvar, muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(S *
    Xvar * ((0.5 * .e22 - (((0.5 * .e18 + .e5)/.e6^2 - .e19 * .e8/((.e3) * .e9 *
    .e5)) * .e1 - (.e19 * .e7/.e2 + 1/.e5)/(.e3)) * .e8/.e9) * .e1) * wHvar,
    uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(S * Xvar * ((0.5 * .e22 - (.e23 * .e1 + .e7 * .e16/((.e3) *
    .e2)) * .e8/.e9) * .e2) * wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(muHvar *
    (.e12 * (.e12/(.e11 * .e13)^2 + mu/(.e11^3 * .e13)) - (((.e14/.e14 - 1) *
      .e10^2/(.e3) + 1) * .e14/.e14 + (.e7/(.e1 * .e9 * .e5) + .e8 * .e2/(.e9 *
      .e5)^2) * .e8 * .e2/(.e3))/(.e3)) * wHvar, muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(muHvar * ((0.5 * .e24 - ((0.5 * .e18 + .e5) * .e2/.e6^2 -
    .e19 * (.e7/.e1 + .e8 * .e2/(.e9 * .e5))/(.e3)) * .e8/.e9) * .e1 + 0.5 *
    (((1 - mu^2/.e11^2)/(.e11 * .e13) - mu * .e12/(.e11 * .e13)^2) * .e12)) *
    wHvar, uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(muHvar * ((((1/.e5 - .e7 * .e16/.e1)/(.e3) -
    .e23 * .e2) * .e8/.e9 + 0.5 * .e24) * .e2) * wHvar, vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(uHvar * (((.e28 * .e1 + 0.5 *
    .e17 - 0.5)/(.e3) - (((0.5 * .e18 + .e5) * (mu * .e2 - (2 * .e29 + 3 * (S *
    epsilon)) * .e1) + (0.5 * (.e1/(.e3)) - 0.5 * (0.5 * (1 - .e1/(.e3)) + .e1/(.e3))) *
    (1 - .e1/(.e3)) * .e2 * .e7/.e5)/.e6^2 + .e19^2 * (.e7/.e6 + .e8/.e9) * .e1 +
    S * epsilon/.e6) * .e8/.e9) * .e1 + 0.5 * (mu * (0.5 * (mu^2/(.e11^3 * .e13)) -
    (0.5 * (.e11 * .e13) - 0.5 * (mu * .e12))/(.e11 * .e13)^2) * .e12)) * wHvar,
    uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(uHvar *
    (((.e19 * (.e7/.e6 + .e8/.e9) * .e16 - ((0.5 * ((1 - .e1/(.e3)) * .e2/(.e3)) +
      0.5 * ((.e1/(.e3) - 1) * .e2/(.e3) + 1 - 0.5 * ((1 - .e1/(.e3)) * (1 -
        .e2/(.e3))))) * .e7/.e5 + mu * (0.5 * .e18 + .e5) - .e15 * (2 * .e29 +
      S * epsilon))/.e6^2) * .e8/.e9 + .e28/(.e3)) * .e1 * .e2) * wHvar, vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(vHvar *
    (((.e28 * .e2 + 0.5 * .e17 - 0.5)/(.e3) + .e8 * (mu/.e6 - ((((3 * (mu) -
      2 * (.e15 * (.e3) * .e7 * .e5/.e6^2)) * .e2 - S * .e1 * epsilon) * .e15 +
      (0.5 * (.e2/(.e3)) - 0.5 * (0.5 * (1 - .e2/(.e3)) + .e2/(.e3))) * (1 -
        .e2/(.e3)) * .e1 * .e7/.e5)/.e6^2 + (.e7/.e6 + .e8/.e9) * .e2 * .e16^2))/.e9) *
      .e2) * wHvar, vHvar)
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
