################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: generalized exponential - normal                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for generalized exponential-normal distribution
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
cgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  ll <- (log(2) - 1/2 * Wu + log(exp(A) * pnorm(a) - exp(B) *
    pnorm(b)))
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for generalized exponential-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
cstgenexponorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  if (S * m3 > 0) {
    varu <- (abs((-S * m3/9)))^(2/3)
  } else {
    varu <- (-S * m3/9)^(2/3)
  }
  if (m2 < varu) {
    varv <- abs(m2 - 5/4 * varu)
  } else {
    varv <- m2 - 5/4 * varu
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
    beta <- c(olsObj[1] + S * sqrt(varu) * 3/2, olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------
#' gradient for generalized exponential-normal distribution
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
cgradgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- exp(Wv/2)
  .e4 <- exp(Wu/2)
  .e5 <- (.e3/.e4 + S * epsilon/.e3)
  .e6 <- pnorm(-.e5)
  .e7 <- dnorm(-.e5)
  .e8 <- (.e3/.e4)
  .e9 <- (2 * .e8 + S * epsilon/.e3)
  .e10 <- dnorm(-.e9)
  .e11 <- pnorm(-.e9)
  .e12 <- exp(.e2/(2 * .e1) + S * epsilon/.e4)
  .e13 <- exp(2 * (.e2/.e1) + 2 * (S * epsilon/.e4))
  .e14 <- (.e10/.e3 - 2 * (.e11/.e4))
  .e15 <- ((.e7/.e3 - .e6/.e4) * .e12 - .e14 * .e13)
  .e16 <- (.e12 * .e6 - .e13 * .e11)
  .e17 <- (.e1 * .e2/(2 * .e1)^2)
  .e18 <- (0.5 * (S * epsilon/.e4) + 2 * .e17)
  .e19 <- (2 * (.e2/.e1) + S * epsilon/.e4)
  .e20 <- (.e10 * .e3/.e4 - .e19 * .e11)
  .e21 <- (0.5 * (.e7 * .e3/.e4) - .e18 * .e6)
  .e22 <- (.e21 * .e12 - .e20 * .e13)
  .e23 <- (0.5 * .e8 - 0.5 * (S * epsilon/.e3))
  .e24 <- (.e2 * .e6/(2 * .e1) - .e23 * .e7)
  .e25 <- (.e3/.e4 - 0.5 * (S * epsilon/.e3))
  .e26 <- (2 * (.e2 * .e11/.e1) - .e10 * .e25)
  .e27 <- (.e12 * .e24 - .e26 * .e13)
  gradll <- cbind(S * Xvar * .e15/.e16, uHvar * (.e22/.e16 - 0.5), vHvar * .e27/.e16)
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for generalized exponential-normal distribution
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
chessgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- exp(Wu)
  .e2 <- exp(Wv)
  .e3 <- exp(Wv/2)
  .e4 <- exp(Wu/2)
  .e5 <- (.e3/.e4 + S * epsilon/.e3)
  .e6 <- pnorm(-.e5)
  .e7 <- dnorm(-.e5)
  .e8 <- (.e3/.e4)
  .e9 <- (2 * .e8 + S * epsilon/.e3)
  .e10 <- dnorm(-.e9)
  .e11 <- pnorm(-.e9)
  .e12 <- exp(.e2/(2 * .e1) + S * epsilon/.e4)
  .e13 <- exp(2 * (.e2/.e1) + 2 * (S * epsilon/.e4))
  .e14 <- (.e10/.e3 - 2 * (.e11/.e4))
  .e15 <- ((.e7/.e3 - .e6/.e4) * .e12 - .e14 * .e13)
  .e16 <- (.e12 * .e6 - .e13 * .e11)
  .e17 <- (.e1 * .e2/(2 * .e1)^2)
  .e18 <- (0.5 * (S * epsilon/.e4) + 2 * .e17)
  .e19 <- (2 * (.e2/.e1) + S * epsilon/.e4)
  .e20 <- (.e10 * .e3/.e4 - .e19 * .e11)
  .e21 <- (0.5 * (.e7 * .e3/.e4) - .e18 * .e6)
  .e22 <- (.e21 * .e12 - .e20 * .e13)
  .e23 <- (0.5 * .e8 - 0.5 * (S * epsilon/.e3))
  .e24 <- (.e2 * .e6/(2 * .e1) - .e23 * .e7)
  .e25 <- (.e3/.e4 - 0.5 * (S * epsilon/.e3))
  .e26 <- (2 * (.e2 * .e11/.e1) - .e10 * .e25)
  .e27 <- (.e12 * .e24 - .e26 * .e13)
  hessll <- matrix(0, nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar + nuZUvar +
    nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * ((((.e5/.e3 - 1/.e4) * .e7/.e3 -
    (.e7/.e3 - .e6/.e4)/.e4) * .e12 - (((.e9/.e3 - 2/.e4) * .e10/.e3 - 2 * (.e14/.e4)) *
    .e13 + .e15^2/.e16))/.e16) * wHvar, Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(S * Xvar * (((((0.5 +
    0.5 * (S * epsilon/.e4) + 2 * .e17) * .e6 - 0.5 * (.e7 * .e3/.e4))/.e4 +
    (0.5 * (.e5/.e4) - .e18/.e3) * .e7) * .e12 - (((.e9/.e4 - .e19/.e3) * .e10 +
    (.e11 - 2 * .e20)/.e4) * .e13 + .e22 * .e15/.e16))/.e16) * wHvar, uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(S *
    Xvar * (((.e7 * (.e2/(2 * .e1) - (.e23 * .e5 + 0.5))/.e3 - .e24/.e4) * .e12 -
    (((2 * (.e2/.e1) - (.e9 * .e25 + 0.5)) * .e10/.e3 - 2 * (.e26/.e4)) * .e13 +
      .e15 * .e27/.e16))/.e16) * wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    ((((0.5 * (0.5 * (.e3 * .e5/.e4) - 0.5) - 0.5 * .e18) * .e7 * .e3/.e4 - (.e21 *
      .e18 + (2 * ((1 - 8 * (.e1^2/(2 * .e1)^2)) * .e1 * .e2/(2 * .e1)^2) -
      0.25 * (S * epsilon/.e4)) * .e6)) * .e12 - ((((.e9 * .e3 - S * epsilon)/.e4 -
      (0.5 + 2 * (.e2/.e1))) * .e10 * .e3/.e4 + (0.5 * (S * epsilon/.e4) +
      2 * (.e2/.e1)) * .e11 - .e19 * .e20) * .e13 + .e22^2/.e16))/.e16) * wHvar,
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * ((((.e7 * .e3/(4 * (.e1 * .e4)) - 2 * (.e1 *
    .e6/(2 * .e1)^2)) * .e2 - ((0.5 * (.e23 * .e5) - 0.25) * .e7 * .e3/.e4 +
    .e18 * .e24)) * .e12 - (.e22 * .e27/.e16 + (2 * ((.e10 * .e3/.e4 - .e11) *
    .e2/.e1) - ((.e9 * .e25 - 0.5) * .e10 * .e3/.e4 + .e26 * .e19)) * .e13))/.e16) *
    wHvar, vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * ((((.e24/2 + (.e6 -
    .e23 * .e7)/2) * .e2/.e1 - (0.25 * .e8 + 0.25 * (S * epsilon/.e3) - .e23^2 *
    .e5) * .e7) * .e12 - (((2 * .e26 + 2 * (.e11 - .e10 * .e25)) * .e2/.e1 -
    (0.25 * (S * epsilon/.e3) + 0.5 * .e8 - .e9 * .e25^2) * .e10) * .e13 + .e27^2/.e16))/.e16) *
    wHvar, vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for generalized exponential-normal distribution
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
genexponormAlgOpt <- function(start, randStart, sdStart, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar, method, derivs, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac, accuracy, stepsize) {
  ## starting values and log likelihood ------
  startVal <- if (!is.null(start))
    start else cstgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], S = S,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar)
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cgenexponormlike(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    cgradbhhhgenexponormAD <- RTMB::MakeADFun(function(p) sum(cgenexponormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = TRUE, silent = TRUE)
    cgradhessgenexponormAD <- RTMB::MakeADFun(function(p) sum(cgenexponormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = FALSE, silent = TRUE)
    fnGradObsgenexpo <- function(parm) {
      Ta1 <- RTMB::GetTape(cgradbhhhgenexponormAD)
      Ta2 <- RTMB::MakeTape(function(weight) {
        WT <- RTMB::MakeTape(function(x) sum(Ta1(x) * weight), parm)
        (WT$jacfun())(RTMB::advector(parm))
      }, rep(1, nrow(Xvar)))
      t((Ta2$jacfun())(RTMB::advector(rep(1, nrow(Xvar)))))
    }
    ### solve for different algorithms ------
    mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(p) -cgradhessgenexponormAD$fn(p),
      gr = function(p) -cgradhessgenexponormAD$gr(p), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
        maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
      maxLikAlgo = maxRoutine(fn = cgradhessgenexponormAD$fn, grad = cgradhessgenexponormAD$gr,
        hess = cgradhessgenexponormAD$he, start = startVal, finalHessian = TRUE,
        control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax,
          reltol = tol, tol = tol, qac = qac)), bhhh = maxLik::maxBHHH(fn = cgradbhhhgenexponormAD$fn,
        grad = fnGradObsgenexpo, hess = cgradhessgenexponormAD$he, start = startVal,
        finalHessian = TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac)), sr1 = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhessgenexponormAD$fn(p), gr = function(p) -cgradhessgenexponormAD$gr(p),
        method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhessgenexponormAD$fn(p), gr = function(p) -cgradhessgenexponormAD$gr(p),
        hs = function(p) as(-cgradhessgenexponormAD$he(p), "dgCMatrix"), method = "Sparse",
        control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(p) -cgradhessgenexponormAD$fn(p),
        gr = function(p) -cgradhessgenexponormAD$gr(p), hess = function(p) -cgradhessgenexponormAD$he(p),
        print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
      nlminb = nlminb(start = startVal, objective = function(p) -cgradhessgenexponormAD$fn(p),
        gradient = function(p) -cgradhessgenexponormAD$gr(p), hessian = function(p) -cgradhessgenexponormAD$he(p),
        control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
          eval.max = itermax, rel.tol = tol, x.tol = tol)))
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$gradient <- cgradhessgenexponormAD$gr(mleObj$par)
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
        mleObj$hessian <- cgradhessgenexponormAD$he(mleObj$par)
      if (method == "sr1")
        mleObj$hessian <- cgradhessgenexponormAD$he(mleObj$solution)
      if (method == "mla")
        mleObj$hessian <- cgradhessgenexponormAD$he(mleObj$b)
    }
    mleObj$logL_OBS <- cgradbhhhgenexponormAD$fn(mlParam)
    mleObj$gradL_OBS <- fnGradObsgenexpo(mlParam)
    rm(cgradbhhhgenexponormAD, cgradhessgenexponormAD, fnGradObsgenexpo)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
    }
    ## numerical derivatives ------
    if (derivs == "numerical") {
      cgradgenexponormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::jacobian(cgenexponormlike, var = unname(parm), params = list(nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar), accuracy = accuracy,
          stepsize = stepsize)
      }
      chessgenexponormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::hessian(function(parm) sum(cgenexponormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)), var = unname(parm),
          accuracy = accuracy, stepsize = stepsize)
      }
      ### solve for different algorithms ------
      mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
        gr = function(parm) -colSums(cgradgenexponormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
          maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
        maxLikAlgo = maxRoutine(fn = cgenexponormlike, grad = cgradgenexponormlikeNum,
          hess = chessgenexponormlikeNum, start = startVal, finalHessian = if (hessianType ==
          2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), method = "SR1", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hs = function(parm) as(-chessgenexponormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"), method = "Sparse",
          control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cgenexponormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          gr = function(parm) -colSums(cgradgenexponormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hess = function(parm) -chessgenexponormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo,
          maxiter = itermax, epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradgenexponormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = function(parm) -chessgenexponormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
      if (method %in% c("ucminf", "nlminb")) {
        mleObj$gradient <- colSums(cgradgenexponormlikeNum(mleObj$par, nXvar = nXvar,
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
          mleObj$hessian <- chessgenexponormlikeNum(parm = mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        if (method == "sr1")
          mleObj$hessian <- chessgenexponormlikeNum(parm = mleObj$solution,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
      mleObj$logL_OBS <- cgenexponormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      mleObj$gradL_OBS <- cgradgenexponormlikeNum(parm = mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
    } else {
      ## analytical derivatives ------
      if (derivs == "analytical") {
        ### solve for different algorithms ------
        mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
          fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
          stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cgenexponormlike,
          grad = cgradgenexponormlike, hess = chessgenexponormlike, start = startVal,
          finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hs = function(parm) as(-chessgenexponormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"),
          method = "Sparse", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
          fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hess = function(parm) -chessgenexponormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo, maxiter = itermax,
          epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradgenexponormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = function(parm) -chessgenexponormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
        if (method %in% c("ucminf", "nlminb")) {
          mleObj$gradient <- colSums(cgradgenexponormlike(mleObj$par, nXvar = nXvar,
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
          mleObj$hessian <- chessgenexponormlike(parm = mleObj$par, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
          if (method == "sr1")
          mleObj$hessian <- chessgenexponormlike(parm = mleObj$solution, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        }
        mleObj$logL_OBS <- cgenexponormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        mleObj$gradL_OBS <- cgradgenexponormlike(parm = mlParam, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
    }
  }
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for generalized exponential-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
cgenexponormeff <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  u <- exp(Wv/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) *
    (dnorm(b) + b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) *
    pnorm(b))
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- (exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a -
      exp(Wv/2)) - exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) *
      pnorm(b - exp(Wv/2)))/(exp(A) * pnorm(a) - exp(B) *
      pnorm(b))
    teBC_reciprocal <- (exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) *
      pnorm(a + exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) +
      exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(exp(A) * pnorm(a) -
      exp(B) * pnorm(b))
    res <- data.frame(u = u, teJLMS = teJLMS, teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for generalized exponential-normal distribution
#' @param object object of class sfacross
#' @noRd
cmarggenexponorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}

cmarggenexponorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
