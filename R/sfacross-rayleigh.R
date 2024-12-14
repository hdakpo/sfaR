################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: rayleigh - normal                                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for rayleigh-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @noRd
craynormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  ll <- (-Wu - 1/2 * Wv - (S * epsilon)^2/(2 * exp(Wv)) + 1/2 *
    (mustar/sigmastar)^2 + 1/2 * log(sigmastar^2) + log(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar)))
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for rayleigh-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
cstraynorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
  vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  if (S * m3 > 0) {
    varu <- exp(1/2 * log(((2 * (S * m3)^2/(pi * (pi - 3)^2))^2)^(2/6)))
  } else {
    varu <- (2 * m3^2/(pi * (3 - pi)^2))^(2/6)
  }
  if (m2 < ((4 - pi)/2) * varu) {
    varv <- abs(2 - ((4 - pi)/2) * varu)
  } else {
    varv <- m2 - ((4 - pi)/2) * varu
  }
  dep_u <- 1/2 * log(((epsiRes^2 - varv) * 2/(4 - pi))^2)
  dep_v <- 1/2 * log((epsiRes^2 - (4 - pi/2) * varu)^2)
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
    beta <- c(olsObj[1] + S * sqrt(varu * pi/2), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------
#' gradient for rayleigh-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @noRd
cgradraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  sigma_sq <- exp(Wu) + exp(Wv)
   .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- (.e1 + .e2)
  .e4 <- sqrt(.e1 * .e2/.e3)
  .e5 <- (.e3 * .e4)
  .e6 <- (S * .e1 * epsilon/.e5)
  .e7 <- pnorm(-.e6)
  .e8 <- dnorm(-.e6)
  .e9 <- (.e8 * .e4 - S * .e1 * .e7 * epsilon/.e3)
  .e10 <- (.e7 + S * .e8 * (.e4/.e2 - .e1/.e5) * epsilon)
  .e11 <- (0.5 * (.e8 * .e2/.e4) - S * .e7 * epsilon)
  .e12 <- (1/.e5 - (0.5 * ((1 - .e1/.e3) * .e2/.e4) + .e4) * .e1/.e5^2)
  .e13 <- (1 - .e2/.e3)
  .e14 <- (0.5 * (.e13 * .e8/.e4) + S * .e7 * epsilon/.e3)
  .e15 <- (0.5 * (.e13 * .e1/.e4) + .e4)
  .e16 <- (.e14/.e9 - .e15 * .e1 * epsilon^2/(.e5^2 * .e4))
  gradll <- cbind(S * Xvar * ((.e10/.e9 - S * epsilon/.e2) * .e1/.e3 + S * epsilon/.e2),
    uHvar * ((.e11 * .e1/(.e9 * .e3) + 0.5) * (1 - .e1/.e3) + .e12 * .e1^2 *
      epsilon^2/.e5 - 1), vHvar * ((.e16 * .e1/.e3 + 2 * ((S * epsilon)^2/(2 *
      .e2)^2)) * .e2 + 0.5 * .e13 - 0.5))
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for rayleigh-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @noRd
chessraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  sigma_sq <- exp(Wu) + exp(Wv)
  .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- (.e1 + .e2)
  .e4 <- sqrt(.e1 * .e2/.e3)
  .e5 <- (.e3 * .e4)
  .e6 <- (S * .e1 * epsilon/.e5)
  .e7 <- pnorm(-.e6)
  .e8 <- dnorm(-.e6)
  .e9 <- (.e8 * .e4 - S * .e1 * .e7 * epsilon/.e3)
  .e10 <- (.e7 + S * .e8 * (.e4/.e2 - .e1/.e5) * epsilon)
  .e11 <- (0.5 * (.e8 * .e2/.e4) - S * .e7 * epsilon)
  .e12 <- (1/.e5 - (0.5 * ((1 - .e1/.e3) * .e2/.e4) + .e4) * .e1/.e5^2)
  .e13 <- (1 - .e2/.e3)
  .e14 <- (0.5 * (.e13 * .e8/.e4) + S * .e7 * epsilon/.e3)
  .e15 <- (0.5 * (.e13 * .e1/.e4) + .e4)
  .e16 <- (.e14/.e9 - .e15 * .e1 * epsilon^2/(.e5^2 * .e4))
  hessll <- matrix(0, nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar + nuZUvar +
    nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * (((((.e1 * epsilon^2/(.e3 * .e2) -
    1) * (.e4/.e2 - .e1/.e5) + .e1/.e5) * .e8 - .e1 * .e10^2/(.e9 * .e3))/.e9 +
    1/.e2) * .e1/.e3 - 1/.e2) * wHvar, Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(S * Xvar * ((((.e7 -
    0.5 * (S * .e8 * .e1 * epsilon/.e5))/(.e9 * .e3) - .e11 * .e1 * .e10/(.e9 *
    .e3)^2) * (1 - .e1/.e3) - 2 * (S * .e12 * .e1 * epsilon/.e5)) * .e1) * wHvar,
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(S *
    Xvar * ((((.e1 * (S * (0.5 * (.e13/.e2) + 1/.e3) * .e8 * epsilon/.e4 - .e14 *
    .e10/.e9) - .e7)/(.e9 * .e3) + 2 * (S * .e15 * .e1 * epsilon/(.e5^2 * .e4))) *
    .e1/.e3 - 4 * (S * epsilon/(2 * .e2)^2)) * .e2) * wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    (((((0.5 * (.e8 * .e2/.e4) + S * (S * .e12 * .e8 * .e1 * epsilon - .e7) *
      epsilon - (.e11 * .e1/.e3 + 0.5 * (0.5 * ((1 - .e1/.e3) * .e8 * .e2/.e4) +
      .e12 * .e8 * .e1 * epsilon^2)))/.e9 - 0.5)/.e3 - (.e11 * (1 - .e1/.e3) +
      .e8 * .e4 - S * .e1 * .e7 * epsilon/.e3) * .e11 * .e1/(.e9 * .e3)^2) *
      (1 - .e1/.e3) + ((2 * .e12 - ((0.5 * (.e1/.e3) + 1 - 0.5 * (0.5 * (1 -
      .e1/.e3) + .e1/.e3)) * (1 - .e1/.e3) * .e2/.e4 + (2 - 2 * ((0.5 * ((1 -
      .e1/.e3) * .e2/.e4) + .e4)^2 * .e1 * .e3/.e5^2)) * .e4) * .e1/.e5^2)/.e5 -
      (0.5 * ((1 - .e1/.e3) * .e2/.e4) + .e4) * .e12 * .e1/.e5^2) * .e1 * epsilon^2) *
      .e1) * wHvar, uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * ((((0.5 * (.e13 * .e8) + 0.5 * ((.e8 * .e2/.e3 -
    .e13 * .e12 * .e8 * .e1 * epsilon^2/.e4) * .e1/.e3 - 0.5 * ((1 - .e1/.e3) *
    .e13 * .e8)))/.e4 + (S * .e7 * epsilon - (.e14 * .e11 * (1 - .e1/.e3)/.e9 +
    S * (.e7/.e3 + S * .e12 * .e8 * epsilon) * epsilon) * .e1)/.e3)/.e9 + 0.5/.e3 -
    (.e16/.e3 + (((0.5 * ((1 - .e1/.e3) * .e2/.e3) + 0.5 * ((.e1/.e3 - 1) * .e2/.e3 +
      1 - 0.5 * ((1 - .e1/.e3) * .e13)) + 0.5 * .e13) * .e1/.e4 + .e4)/(.e5^2 *
      .e4) + .e15 * (1/(.e5^2 * .e4) - (0.5 * (.e5^2 * (1 - .e1/.e3)/.e5) +
      2 * ((0.5 * ((1 - .e1/.e3) * .e2/.e4) + .e4) * .e1)) * .e1 * .e2/(.e5^2 *
      .e4)^2)) * epsilon^2) * .e1) * .e1 * .e2/.e3) * wHvar, vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * ((((((0.5 * (.e2 * (.e15 *
    .e8 * .e1^2 * epsilon^2/(.e5^2 * .e4) - .e8)/.e3 - 0.5 * (.e13 * .e8)) +
    0.5 * .e8) * .e13/.e4 + (.e2 * (S * (S * .e15 * .e8 * .e1 * epsilon/.e5^2 -
    .e7/.e3) * epsilon - .e14^2 * .e1/.e9) + S * .e7 * epsilon)/.e3)/.e9 - (.e16 *
    .e2/.e3 + (.e15 * (1/(.e5^2 * .e4) - (0.5 * (.e5^2 * .e13/.e5) + 2 * (.e15 *
    .e2)) * .e1 * .e2/(.e5^2 * .e4)^2) + (0.5 * (.e2/.e3) - 0.5 * (0.5 * .e13 +
    .e2/.e3)) * .e13 * .e3/(.e5^2 * .e2)) * .e1 * epsilon^2)) * .e1 - 0.5 * .e13)/.e3 +
    (2 - 16 * (.e2^2/(2 * .e2)^2)) * (S * epsilon)^2/(2 * .e2)^2) * .e2) * wHvar,
    vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for rayleigh-normal distribution
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
raynormAlgOpt <- function(start, randStart, sdStart, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar, method, derivs, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac, accuracy, stepsize) {
  ## starting values and log likelihood ------
  startVal <- if (!is.null(start))
    start else cstraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], S = S,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar)
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(craynormlike(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    cgradbhhhraynormAD <- RTMB::MakeADFun(function(p) sum(craynormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = TRUE, silent = TRUE)
    cgradhessraynormAD <- RTMB::MakeADFun(function(p) sum(craynormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = FALSE, silent = TRUE)
    fnGradObsray <- function(parm) {
      Ta1 <- RTMB::GetTape(cgradbhhhraynormAD)
      Ta2 <- RTMB::MakeTape(function(weight) {
        WT <- RTMB::MakeTape(function(x) sum(Ta1(x) * weight), parm)
        (WT$jacfun())(RTMB::advector(parm))
      }, rep(1, nrow(Xvar)))
      t((Ta2$jacfun())(RTMB::advector(rep(1, nrow(Xvar)))))
    }
    ### solve for different algorithms ------
    mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(p) -cgradhessraynormAD$fn(p),
      gr = function(p) -cgradhessraynormAD$gr(p), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
        maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
      maxLikAlgo = maxRoutine(fn = cgradhessraynormAD$fn, grad = cgradhessraynormAD$gr,
        hess = cgradhessraynormAD$he, start = startVal, finalHessian = TRUE,
        control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax,
          reltol = tol, tol = tol, qac = qac)), bhhh = maxLik::maxBHHH(fn = cgradbhhhraynormAD$fn,
        grad = fnGradObsray, hess = cgradhessraynormAD$he, start = startVal,
        finalHessian = TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac)), sr1 = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhessraynormAD$fn(p), gr = function(p) -cgradhessraynormAD$gr(p),
        method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhessraynormAD$fn(p), gr = function(p) -cgradhessraynormAD$gr(p),
        hs = function(p) as(-cgradhessraynormAD$he(p), "dgCMatrix"), method = "Sparse",
        control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(p) -cgradhessraynormAD$fn(p),
        gr = function(p) -cgradhessraynormAD$gr(p), hess = function(p) -cgradhessraynormAD$he(p),
        print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
      nlminb = nlminb(start = startVal, objective = function(p) -cgradhessraynormAD$fn(p),
        gradient = function(p) -cgradhessraynormAD$gr(p), hessian = function(p) -cgradhessraynormAD$he(p),
        control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
          eval.max = itermax, rel.tol = tol, x.tol = tol)))
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$gradient <- cgradhessraynormAD$gr(mleObj$par)
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
        mleObj$hessian <- cgradhessraynormAD$he(mleObj$par)
      if (method == "sr1")
        mleObj$hessian <- cgradhessraynormAD$he(mleObj$solution)
      if (method == "mla")
        mleObj$hessian <- cgradhessraynormAD$he(mleObj$b)
    }
    mleObj$logL_OBS <- cgradbhhhraynormAD$fn(mlParam)
    mleObj$gradL_OBS <- fnGradObsray(mlParam)
    rm(cgradbhhhraynormAD, cgradhessraynormAD, fnGradObsray)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
    }
    ## numerical derivatives ------
    if (derivs == "numerical") {
      cgradraynormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::jacobian(craynormlike, var = unname(parm), params = list(nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar), accuracy = accuracy,
          stepsize = stepsize)
      }
      chessraynormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::hessian(function(parm) sum(craynormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)), var = unname(parm),
          accuracy = accuracy, stepsize = stepsize)
      }
      ### solve for different algorithms ------
      mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(craynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
        gr = function(parm) -colSums(cgradraynormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
          maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
        maxLikAlgo = maxRoutine(fn = craynormlike, grad = cgradraynormlikeNum,
          hess = chessraynormlikeNum, start = startVal, finalHessian = if (hessianType ==
          2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(craynormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradraynormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), method = "SR1", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(craynormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradraynormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hs = function(parm) as(-chessraynormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"), method = "Sparse",
          control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(craynormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          gr = function(parm) -colSums(cgradraynormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hess = function(parm) -chessraynormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo,
          maxiter = itermax, epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(craynormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradraynormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = function(parm) -chessraynormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
      if (method %in% c("ucminf", "nlminb")) {
        mleObj$gradient <- colSums(cgradraynormlikeNum(mleObj$par, nXvar = nXvar,
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
          mleObj$hessian <- chessraynormlikeNum(parm = mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        if (method == "sr1")
          mleObj$hessian <- chessraynormlikeNum(parm = mleObj$solution,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
      mleObj$logL_OBS <- craynormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      mleObj$gradL_OBS <- cgradraynormlikeNum(parm = mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
    } else {
      ## analytical derivatives ------
      if (derivs == "analytical") {
        ### solve for different algorithms ------
        mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
          fn = function(parm) -sum(craynormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradraynormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
          stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = craynormlike,
          grad = cgradraynormlike, hess = chessraynormlike, start = startVal,
          finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(craynormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradraynormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(craynormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradraynormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hs = function(parm) as(-chessraynormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"),
          method = "Sparse", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
          fn = function(parm) -sum(craynormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradraynormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hess = function(parm) -chessraynormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo, maxiter = itermax,
          epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(craynormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradraynormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = function(parm) -chessraynormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
        if (method %in% c("ucminf", "nlminb")) {
          mleObj$gradient <- colSums(cgradraynormlike(mleObj$par, nXvar = nXvar,
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
          mleObj$hessian <- chessraynormlike(parm = mleObj$par, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
          if (method == "sr1")
          mleObj$hessian <- chessraynormlike(parm = mleObj$solution, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        }
        mleObj$logL_OBS <- craynormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        mleObj$gradL_OBS <- cgradraynormlike(parm = mlParam, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
    }
  }
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for rayleigh-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
craynormeff <- function(object, level) {
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
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta), t(Xvar))[1, ]
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  u <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 +
    sigmastar^2) * pnorm(mustar/sigmastar))/(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  m <- ifelse(mustar/2 + sqrt(sigmastar^2 + mustar^2/4) > 0,
    mustar/2 + sqrt(sigmastar^2 + mustar^2/4), 0)
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar))/(sigmastar * dnorm(mustar/sigmastar) +
      mustar * pnorm(mustar/sigmastar))
    teBC_reciprocal <- exp(mustar + sigmastar^2/2) * (sigmastar *
      dnorm(mustar/sigmastar + sigmastar) + (mustar + sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teMO <- exp(-m)
    res <- data.frame(u = u, teJLMS = teJLMS, teBC = teBC,
      m = m, teMO = teMO, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u, m = m)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for rayleigh-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargraynorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2) * 1/2 *
    sqrt(pi/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}

cmargraynorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) * (4 -
    pi)/2, ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
