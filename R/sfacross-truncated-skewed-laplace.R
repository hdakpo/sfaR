################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: truncated skewed laplace - normal                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for truncated skewed laplace-normal distribution
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
#' @noRd
ctslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  A <- exp(Wv)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  if (lambda < 0) {
    return(NA)
  } else {
    ll <- log(1 + lambda) - log(2 * lambda + 1) - 1/2 * Wu +
      log(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
    RTMB::ADREPORT(ll * wHvar)
    return(ll * wHvar)
  }
}

# starting value for the log-likelihood ----------
#' starting values for truncated skewed laplace-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
csttslnorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
  vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  if (S * m3 > 0) {
    varu <- (abs((-S * m3/2)))^(2/3)
  } else {
    varu <- (-S * m3/2)^(2/3)
  }
  if (m2 < varu) {
    varv <- abs(m2 - varu)
  } else {
    varv <- m2 - varu
  }
  dep_u <- 1/2 * log((epsiRes^2 - varv)^2)
  dep_v <- 1/2 * log((epsiRes^2 - varu)^2)
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
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi, lambda = 1))
}

# Gradient of the likelihood function ----------
#' gradient for truncated skewed laplace-normal distribution
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
#' @noRd
cgradtslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- exp(Wu/2)
  .e4 <- exp(Wv/2)
  .e5 <- (.e4/.e3 + S * epsilon/.e4)
  .e6 <- pnorm(-.e5)
  .e7 <- dnorm(-.e5)
  .e8 <- ((1 + lambda) * .e4/.e3 + S * epsilon/.e4)
  .e9 <- pnorm(-.e8)
  .e10 <- dnorm(-.e8)
  .e11 <- exp(.e2/(2 * .e1) + S * epsilon/.e3)
  .e12 <- exp(((1 + lambda) * .e2/(2 * .e1) + S * epsilon/.e3) * (1 + lambda))
  .e13 <- (.e10/.e4 - (1 + lambda) * .e9/.e3)
  .e14 <- (.e7/.e4 - .e6/.e3)
  .e15 <- (2 * (.e14 * .e11) - .e13 * .e12)
  .e16 <- (2 * (.e11 * .e6) - .e12 * .e9)
  .e17 <- (.e1 * .e2/(2 * .e1)^2)
  .e18 <- (0.5 * (S * epsilon/.e3) + 2 * .e17)
  .e19 <- ((1 + lambda) * .e1 * .e2/(2 * .e1)^2)
  .e20 <- (0.5 * (S * epsilon/.e3) + 2 * .e19)
  .e21 <- (0.5 * (.e10 * .e4/.e3) - .e20 * .e9)
  .e22 <- ((0.5 * (.e7 * .e4/.e3) - .e18 * .e6) * .e11)
  .e23 <- (2 * .e22 - .e21 * (1 + lambda) * .e12)
  .e24 <- (0.5 * (.e4/.e3) - 0.5 * (S * epsilon/.e4))
  .e25 <- (.e11 * (.e2 * .e6/(2 * .e1) - .e24 * .e7))
  .e26 <- (0.5 * ((1 + lambda) * .e4/.e3) - 0.5 * (S * epsilon/.e4))
  .e27 <- (2 * .e25 - ((1 + lambda)^2 * .e2 * .e9/(2 * .e1) - .e26 * .e10) * .e12)
  .e28 <- ((1 + lambda) * .e2/.e1 + S * epsilon/.e3)
  .e29 <- (.e28 * .e9 - .e10 * .e4/.e3)
  gradll <- cbind(S * Xvar * .e15/.e16, uHvar * (.e23/.e16 - 0.5), vHvar * .e27/.e16,
    1/(1 + lambda) - (.e29 * .e12/.e16 + 2/(1 + 2 * lambda)))
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for truncated skewed laplace-normal distribution
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
#' @noRd
chesstslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
 .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- exp(Wu/2)
  .e4 <- exp(Wv/2)
  .e5 <- (.e4/.e3 + S * epsilon/.e4)
  .e6 <- pnorm(-.e5)
  .e7 <- dnorm(-.e5)
  .e8 <- ((1 + lambda) * .e4/.e3 + S * epsilon/.e4)
  .e9 <- pnorm(-.e8)
  .e10 <- dnorm(-.e8)
  .e11 <- exp(.e2/(2 * .e1) + S * epsilon/.e3)
  .e12 <- exp(((1 + lambda) * .e2/(2 * .e1) + S * epsilon/.e3) * (1 + lambda))
  .e13 <- (.e10/.e4 - (1 + lambda) * .e9/.e3)
  .e14 <- (.e7/.e4 - .e6/.e3)
  .e15 <- (2 * (.e14 * .e11) - .e13 * .e12)
  .e16 <- (2 * (.e11 * .e6) - .e12 * .e9)
  .e17 <- (.e1 * .e2/(2 * .e1)^2)
  .e18 <- (0.5 * (S * epsilon/.e3) + 2 * .e17)
  .e19 <- ((1 + lambda) * .e1 * .e2/(2 * .e1)^2)
  .e20 <- (0.5 * (S * epsilon/.e3) + 2 * .e19)
  .e21 <- (0.5 * (.e10 * .e4/.e3) - .e20 * .e9)
  .e22 <- ((0.5 * (.e7 * .e4/.e3) - .e18 * .e6) * .e11)
  .e23 <- (2 * .e22 - .e21 * (1 + lambda) * .e12)
  .e24 <- (0.5 * (.e4/.e3) - 0.5 * (S * epsilon/.e4))
  .e25 <- (.e11 * (.e2 * .e6/(2 * .e1) - .e24 * .e7))
  .e26 <- (0.5 * ((1 + lambda) * .e4/.e3) - 0.5 * (S * epsilon/.e4))
  .e27 <- (2 * .e25 - ((1 + lambda)^2 * .e2 * .e9/(2 * .e1) - .e26 * .e10) * .e12)
  .e28 <- ((1 + lambda) * .e2/.e1 + S * epsilon/.e3)
  .e29 <- (.e28 * .e9 - .e10 * .e4/.e3)
  hessll <- matrix(0, nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar + nuZUvar +
    nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * ((2 * (((.e5/.e4 - 1/.e3) * .e7/.e4 -
    .e14/.e3) * .e11) - (((.e8/.e4 - (1 + lambda)/.e3) * .e10/.e4 - (1 + lambda) *
    .e13/.e3) * .e12 + .e15^2/.e16))/.e16) * wHvar, Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(S * Xvar * ((2 *
    ((((0.5 + 0.5 * (S * epsilon/.e3) + 2 * .e17) * .e6 - 0.5 * (.e7 * .e4/.e3))/.e3 +
      (0.5 * (.e5/.e3) - .e18/.e4) * .e7) * .e11) - (((0.5 * (.e8/.e3) - .e20/.e4) *
    .e10 + (0.5 * .e9 - .e21 * (1 + lambda))/.e3) * (1 + lambda) * .e12 + .e23 *
    .e15/.e16))/.e16) * wHvar, uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(S *
    Xvar * ((2 * ((.e7 * (.e2/(2 * .e1) - (.e24 * .e5 + 0.5))/.e4 - (.e2 * .e6/(2 *
    .e1) - .e24 * .e7)/.e3) * .e11) - ((((1 + lambda)^2 * .e2/(2 * .e1) - (.e8 *
    .e26 + 0.5)) * .e10/.e4 - ((1 + lambda)^2 * .e2 * .e9/(2 * .e1) - .e26 *
    .e10) * (1 + lambda)/.e3) * .e12 + .e15 * .e27/.e16))/.e16) * wHvar, vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- crossprod(S * Xvar, (-(((.e28/.e4 -
    .e8/.e3) * .e10 - ((.e29 * (1 + lambda) + .e9)/.e3 + .e29 * .e15/.e16)) *
    .e12/.e16)) * wHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    ((2 * (((0.5 * (0.5 * (.e4 * .e5/.e3) - 0.5) - 0.5 * .e18) * .e7 * .e4/.e3 -
      ((0.5 * (.e7 * .e4/.e3) - .e18 * .e6) * .e18 + (2 * ((1 - 8 * (.e1^2/(2 *
        .e1)^2)) * .e1 * .e2/(2 * .e1)^2) - 0.25 * (S * epsilon/.e3)) * .e6)) *
      .e11) - (((0.5 * (0.5 * (.e8 * (1 + lambda) * .e4/.e3) - 0.5) - 0.5 *
      (.e20 * (1 + lambda))) * .e10 * .e4/.e3 - (.e21 * .e20 * (1 + lambda) +
      (2 * ((1 - 8 * (.e1^2/(2 * .e1)^2)) * (1 + lambda) * .e1 * .e2/(2 * .e1)^2) -
        0.25 * (S * epsilon/.e3)) * .e9)) * (1 + lambda) * .e12 + .e23^2/.e16))/.e16) *
    wHvar, uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * ((2 * (((.e7 * .e4/(4 * (.e1 * .e3)) - 2 *
    (.e1 * .e6/(2 * .e1)^2)) * .e2 - ((0.5 * (.e24 * .e5) - 0.25) * .e7 * .e4/.e3 +
    .e18 * (.e2 * .e6/(2 * .e1) - .e24 * .e7))) * .e11) - ((((1 + lambda) * .e10 *
    .e4/(4 * (.e1 * .e3)) - 2 * (.e1 * .e9/(2 * .e1)^2)) * (1 + lambda) * .e2 -
    (((1 + lambda)^2 * .e2 * .e9/(2 * .e1) - .e26 * .e10) * .e20 + (0.5 * (.e8 *
      .e26) - 0.25) * .e10 * .e4/.e3)) * (1 + lambda) * .e12 + .e23 * .e27/.e16))/.e16) *
    wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- crossprod(uHvar,
    (-((((0.5 * .e28 - 0.5 * (.e8 * .e4/.e3)) * (1 + lambda) + 0.5) * .e10 *
      .e4/.e3 - (.e29 * (.e20 * (1 + lambda) + .e23/.e16) + ((1 + lambda) *
      .e2/.e1 + 0.5 * (S * epsilon/.e3)) * .e9)) * .e12/.e16)) * wHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * ((2 * ((((.e2 * .e6/(2 *
    .e1) - .e24 * .e7)/2 + (.e6 - .e24 * .e7)/2) * .e2/.e1 - (0.25 * (.e4/.e3) +
    0.25 * (S * epsilon/.e4) - .e24^2 * .e5) * .e7) * .e11) - (((((1 + lambda)^2 *
    .e2 * .e9/(2 * .e1) - .e26 * .e10)/2 + (.e9 - .e26 * .e10)/2) * (1 + lambda)^2 *
    .e2/.e1 - (0.25 * ((1 + lambda) * .e4/.e3) + 0.25 * (S * epsilon/.e4) - .e8 *
    .e26^2) * .e10) * .e12 + .e27^2/.e16))/.e16) * wHvar, vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- crossprod(vHvar, (-(((.e29 * (1 + lambda)/2 + .e9) * (1 +
    lambda) * .e2/.e1 - ((.e28 * .e26 + (0.5 - .e8 * .e26) * .e4/.e3) * .e10 +
    .e29 * .e27/.e16)) * .e12/.e16)) * wHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum((lambda = 4/(1 +
    2 * lambda)^2 - (((.e29 * .e12/.e16 + (1 + lambda) * .e2/.e1 + S * epsilon/.e3) *
    .e29 + ((.e8 * .e4 - S * epsilon)/.e3 - (1 + lambda) * .e2/.e1) * .e10 *
    .e4/.e3 + .e2 * .e9/.e1) * .e12/.e16 + 1/(1 + lambda)^2)) * wHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for truncated skewed laplace-normal distribution
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
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
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
tslnormAlgOpt <- function(start, randStart, sdStart, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar, method, derivs, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac, accuracy, stepsize) {
  ## starting values and log likelihood ------
  startVal <- if (!is.null(start))
    start else csttslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], S = S,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar)
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ctslnormlike(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S))
  ## automatic differentiation ------
  if (derivs == "ad") {
    if (method %in% c("bfgs", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
    cgradbhhhtslnormAD <- RTMB::MakeADFun(function(p) sum(ctslnormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = TRUE, silent = TRUE)
    cgradhesstslnormAD <- RTMB::MakeADFun(function(p) sum(ctslnormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = FALSE, silent = TRUE)
    fnGradObstsl <- function(parm) {
      Ta1 <- RTMB::GetTape(cgradbhhhtslnormAD)
      Ta2 <- RTMB::MakeTape(function(weight) {
        WT <- RTMB::MakeTape(function(x) sum(Ta1(x) * weight), parm)
        (WT$jacfun())(RTMB::advector(parm))
      }, rep(1, nrow(Xvar)))
      t((Ta2$jacfun())(RTMB::advector(rep(1, nrow(Xvar)))))
    }
    ### solve for different algorithms ------
    mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(p) -cgradhesstslnormAD$fn(p),
      gr = function(p) -cgradhesstslnormAD$gr(p), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
        maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
      maxLikAlgo = maxRoutine(fn = cgradhesstslnormAD$fn, grad = cgradhesstslnormAD$gr,
        hess = cgradhesstslnormAD$he, start = startVal, finalHessian = TRUE,
        control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax,
          reltol = tol, tol = tol, qac = qac)), bhhh = maxLik::maxBHHH(fn = cgradbhhhtslnormAD$fn,
        grad = fnGradObstsl, hess = cgradhesstslnormAD$he, start = startVal,
        finalHessian = TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac)), sr1 = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesstslnormAD$fn(p), gr = function(p) -cgradhesstslnormAD$gr(p),
        method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesstslnormAD$fn(p), gr = function(p) -cgradhesstslnormAD$gr(p),
        hs = function(p) as(-cgradhesstslnormAD$he(p), "dgCMatrix"), method = "Sparse",
        control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(p) -cgradhesstslnormAD$fn(p),
        gr = function(p) -cgradhesstslnormAD$gr(p), hess = function(p) -cgradhesstslnormAD$he(p),
        print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
      nlminb = nlminb(start = startVal, objective = function(p) -cgradhesstslnormAD$fn(p),
        gradient = function(p) -cgradhesstslnormAD$gr(p), hessian = function(p) -cgradhesstslnormAD$he(p),
        control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
          eval.max = itermax, rel.tol = tol, x.tol = tol)))
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$gradient <- cgradhesstslnormAD$gr(mleObj$par)
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
        mleObj$hessian <- cgradhesstslnormAD$he(mleObj$par)
      if (method == "sr1")
        mleObj$hessian <- cgradhesstslnormAD$he(mleObj$solution)
      if (method == "mla")
        mleObj$hessian <- cgradhesstslnormAD$he(mleObj$b)
    }
    mleObj$logL_OBS <- cgradbhhhtslnormAD$fn(mlParam)
    mleObj$gradL_OBS <- fnGradObstsl(mlParam)
    rm(cgradbhhhtslnormAD, cgradhesstslnormAD, fnGradObstsl)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
    }
    ## numerical derivatives ------
    if (derivs == "numerical") {
      cgradtslnormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::jacobian(ctslnormlike, var = unname(parm), params = list(nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar), accuracy = accuracy,
          stepsize = stepsize)
      }
      chesstslnormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::hessian(function(parm) sum(ctslnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)), var = unname(parm),
          accuracy = accuracy, stepsize = stepsize)
      }
      ### solve for different algorithms ------
      mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ctslnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
        gr = function(parm) -colSums(cgradtslnormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
          maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
        maxLikAlgo = maxRoutine(fn = ctslnormlike, grad = cgradtslnormlikeNum,
          hess = chesstslnormlikeNum, start = startVal, finalHessian = if (hessianType ==
          2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtslnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), method = "SR1", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtslnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hs = function(parm) as(-chesstslnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"), method = "Sparse",
          control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(ctslnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          gr = function(parm) -colSums(cgradtslnormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hess = function(parm) -chesstslnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo,
          maxiter = itermax, epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradtslnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = function(parm) -chesstslnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
      if (method %in% c("ucminf", "nlminb")) {
        mleObj$gradient <- colSums(cgradtslnormlikeNum(mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S))
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
          mleObj$hessian <- chesstslnormlikeNum(parm = mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        if (method == "sr1")
          mleObj$hessian <- chesstslnormlikeNum(parm = mleObj$solution,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
      mleObj$logL_OBS <- ctslnormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      mleObj$gradL_OBS <- cgradtslnormlikeNum(parm = mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
    } else {
      ## analytical derivatives ------
      if (derivs == "analytical") {
        ### solve for different algorithms ------
        mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
          fn = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtslnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
          stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ctslnormlike,
          grad = cgradtslnormlike, hess = chesstslnormlike, start = startVal,
          finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtslnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtslnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hs = function(parm) as(-chesstslnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"),
          method = "Sparse", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
          fn = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtslnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hess = function(parm) -chesstslnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo, maxiter = itermax,
          epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradtslnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = function(parm) -chesstslnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
        if (method %in% c("ucminf", "nlminb")) {
          mleObj$gradient <- colSums(cgradtslnormlike(mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S))
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
          mleObj$hessian <- chesstslnormlike(parm = mleObj$par, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
          if (method == "sr1")
          mleObj$hessian <- chesstslnormlike(parm = mleObj$solution, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        }
        mleObj$logL_OBS <- ctslnormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        mleObj$gradL_OBS <- cgradtslnormlike(parm = mlParam, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
    }
  }
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for truncated skewed laplace-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
ctslnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    object$nvZVvar + 1]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta), t(Xvar))[1, ]
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  u <- exp(Wv/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) -
    exp(B) * (dnorm(b) + b * pnorm(b)))/(2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- (2 * exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) *
      pnorm(a - exp(Wv/2)) - exp(B) * exp(-b * exp(Wv/2) +
      exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(2 * exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal <- (2 * exp(A) * exp(a * exp(Wv/2) +
      exp(Wv)/2) * pnorm(a + exp(Wv/2)) - exp(B) * exp(b *
      exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(2 *
      exp(A) * pnorm(a) - exp(B) * pnorm(b))
    res <- data.frame(u = u, teJLMS = teJLMS, teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for truncated skewed laplace-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargtslnorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    object$nvZVvar + 1]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}

cmargtslnorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    object$nvZVvar + 1]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
