################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for halfnormal-normal distribution
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
chalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  ll <- (-log(1/2) - 1/2 * log(exp(Wu) + exp(Wv)) + dnorm(epsilon/sqrt(exp(Wu) +
    exp(Wv)), log = TRUE) + log(pnorm(mustar/sigmastar)))
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for halfnormal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
csthalfnorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
  vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  ## Coelli (1995) suggests 0.05 for gamma when wrong sign
  if (S * m3 > 0) {
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
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu * 2/pi), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------
#' gradient for halfnormal-normal distribution
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
cgradhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  sigma_sq <- exp(Wu) + exp(Wv)
  .e2 <- exp(Wu)
  .e3 <- exp(Wv)
  .e4 <- .e2 + .e3
  .e6 <- epsilon
  .e9 <- sqrt(.e2 * .e3/.e4)
  .e10 <- .e4 * .e9
  .e11 <- -(S * .e2 * .e6/.e10)
  .e13 <- S * .e6/sqrt(.e4)
  .e14 <- dnorm(.e11)
  .e15 <- dnorm(.e13)
  .e17 <- pnorm(.e11)
  .e19 <- S * dnorm(.e13) * .e6
  .e20 <- .e10^2
  .e22 <- 0.5 * (.e19/(.e15 * .e4^2))
  .e23 <- 0.5/.e4
  gradll <- cbind(Xvar * (S * (.e14 * .e2/(.e17 * .e9) + .e19/.e15)/.e4), uHvar *
    (.e2 * (S * (.e22 - (1/.e10 - (0.5 * ((1 - .e2/.e4) * .e3/.e9) + .e9) * .e2/.e20) *
      .e14/.e17) * .e6 - .e23)), vHvar * (.e3 * (S * ((0.5 * ((1 - .e3/.e4) *
    .e2/.e9) + .e9) * .e14 * .e2/(.e20 * .e17) + .e22) * .e6 - .e23)))
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for halfnormal-normal distribution
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
chesshalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  Wv <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  sigma_sq <- exp(Wu) + exp(Wv)
  .e2 <- exp(Wu)
  .e3 <- exp(Wv)
  .e4 <- .e2 + .e3
  .e7 <- sqrt(.e2 * .e3/.e4)
  .e9 <- epsilon
  .e10 <- .e4 * .e7
  .e13 <- S * .e2 * .e9/.e10
  .e14 <- -.e13
  .e15 <- S * .e9
  .e16 <- .e15/sqrt(.e4)
  .e17 <- pnorm(.e14)
  .e18 <- .e10^2
  .e19 <- dnorm(.e16)
  .e20 <- dnorm(.e14)
  .e21 <- .e2/.e4
  .e22 <- .e3/.e4
  .e23 <- S^2
  .e24 <- .e9^2
  .e26 <- 1 - .e21
  .e27 <- 1 - .e22
  .e28 <- .e4^2
  .e29 <- .e26 * .e3
  .e32 <- 0.5 * (.e29/.e7) + .e7
  .e33 <- .e19 * .e28
  .e37 <- 0.5 * (.e27 * .e2/.e7) + .e7
  .e38 <- .e18 * .e17
  .e42 <- 1/.e10 - .e32 * .e2/.e18
  .e43 <- .e17 * .e7
  .e44 <- .e33^2
  .e45 <- .e23 * .e24
  .e46 <- .e23 * .e19
  .e50 <- 0.5 * (.e45/(.e19 * .e4^4)) - (0.5 * (.e46 * .e24) + 2 * (.e19 * .e4))/.e44
  .e51 <- .e2^2
  .e53 <- .e38^2
  .e54 <- .e42 * .e20
  .e55 <- .e43^2
  .e56 <- 0.5 * (S * .e50 * .e19 * .e9)
  .e57 <- 0.5/.e33
  .e58 <- 0.5/.e28
  .e59 <- .e19/.e19
  .e60 <- S * .e19
  .e61 <- .e32 * .e4
  .e62 <- .e37 * .e20
  .e63 <- (0.5 * ((.e45/.e4 - 1)/.e33 - .e46 * .e4 * .e24/.e44) - .e57) * .e19
  .e65 <- .e54/.e17
  .e66 <- .e4 * .e17
  .e68 <- 0.5 * (.e26 * .e27)
  .e69 <- 0.5 * .e27
  .e70 <- 0.5 * .e21
  .e71 <- 0.5 * .e22
  .e74 <- .e20 * .e51
  .e75 <- .e20/.e17
  .e76 <- S * .e42
  .e80 <- .e60 * (.e23 * (0.5 - 0.5 * .e59) * .e24/.e4 - 1) * .e9/.e19 - .e20 *
    .e2/.e43
  .e83 <- .e23 * .e37
  hessll <- matrix(0, nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar + nuZUvar +
    nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * (.e23 * (.e74 * (.e15/(.e3 * .e17 *
    .e7) - .e20/.e55)/.e4 + .e19 * (.e23 * (1 - .e59) * .e24/.e4 - 1)/.e19)/.e4) *
    wHvar, Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(S * Xvar * ((.e65 +
    S * (.e63 - .e54 * .e2 * (.e15/.e3 - .e20/.e43)/.e66) * .e9) * .e2) * wHvar,
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(S *
    Xvar * (.e3 * (S * (.e62 * .e51 * (.e15/(.e18 * .e3 * .e17) - .e18 * .e20/(.e53 *
    .e7))/.e4 + .e63) * .e9 - .e62 * .e2/.e38)) * wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    (((.e58 + S * (.e56 - .e20 * (S * .e42^2 * (.e75 - .e13) * .e9 - ((.e70 +
      1 - 0.5 * (0.5 * .e26 + .e21)) * .e26 * .e3/.e7 + (2 - 2 * (.e32^2 *
      .e2 * .e4/.e18)) * .e7)/.e18)/.e17) * .e9) * .e2 + S * (0.5 * (.e60 *
      .e9/.e33) - .e65) * .e9 - 0.5/.e4) * .e2) * wHvar, uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * ((.e58 + S * (((((0.5 * .e29 - .e83 * .e42 *
    .e2 * .e24)/.e4 + 0.5 * ((.e21 - 1) * .e3/.e4 + 1 - .e68) + .e69) * .e2/.e7 +
    .e7)/.e38 - .e37 * (2 * (.e61 * .e17 * .e7) - S * .e18 * .e42 * .e20 * .e9) *
    .e2/.e53) * .e20 + .e56) * .e9) * .e2 * .e3) * wHvar, vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * (((.e71 - 0.5)/.e4 +
    S * ((((.e71 - 0.5 * (.e69 + .e22)) * .e27 + .e23 * .e37^2 * .e2 * .e3 *
      .e24/(.e18 * .e4)) * .e2/(.e38 * .e7) + .e37 * (1/.e38 - .e37 * (2 *
      (.e66 * .e7) + S * .e20 * .e2 * .e9) * .e3/.e53)) * .e20 * .e2 + S *
      (0.5 * (.e50 * .e3) + .e57) * .e19 * .e9) * .e9) * .e3) * wHvar, vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for halfnormal-normal distribution
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
halfnormAlgOpt <- function(start, randStart, sdStart, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar, method, derivs, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac, accuracy, stepsize) {
  ## starting values and log likelihood ------
  startVal <- if (!is.null(start))
    start else csthalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], S = S,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar)
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(chalfnormlike(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    cgradbhhhhalfnormAD <- RTMB::MakeADFun(function(p) sum(chalfnormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = TRUE, silent = TRUE)
    cgradhesshalfnormAD <- RTMB::MakeADFun(function(p) sum(chalfnormlike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = FALSE, silent = TRUE)
    fnGradObshalf <- function(parm) {
      Ta1 <- RTMB::GetTape(cgradbhhhhalfnormAD)
      Ta2 <- RTMB::MakeTape(function(weight) {
        WT <- RTMB::MakeTape(function(x) sum(Ta1(x) * weight), parm)
        (WT$jacfun())(RTMB::advector(parm))
      }, rep(1, nrow(Xvar)))
      t((Ta2$jacfun())(RTMB::advector(rep(1, nrow(Xvar)))))
    }
    ### solve for different algorithms ------
    mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(p) -cgradhesshalfnormAD$fn(p),
      gr = function(p) -cgradhesshalfnormAD$gr(p), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
        maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
      maxLikAlgo = maxRoutine(fn = cgradhesshalfnormAD$fn, grad = cgradhesshalfnormAD$gr,
        hess = cgradhesshalfnormAD$he, start = startVal, finalHessian = TRUE,
        control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax,
          reltol = tol, tol = tol, qac = qac)), bhhh = maxLik::maxBHHH(fn = cgradbhhhhalfnormAD$fn,
        grad = fnGradObshalf, hess = cgradhesshalfnormAD$he, start = startVal,
        finalHessian = TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac)), sr1 = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesshalfnormAD$fn(p), gr = function(p) -cgradhesshalfnormAD$gr(p),
        method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesshalfnormAD$fn(p), gr = function(p) -cgradhesshalfnormAD$gr(p),
        hs = function(p) as(-cgradhesshalfnormAD$he(p), "dgCMatrix"), method = "Sparse",
        control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(p) -cgradhesshalfnormAD$fn(p),
        gr = function(p) -cgradhesshalfnormAD$gr(p), hess = function(p) -cgradhesshalfnormAD$he(p),
        print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
      nlminb = nlminb(start = startVal, objective = function(p) -cgradhesshalfnormAD$fn(p),
        gradient = function(p) -cgradhesshalfnormAD$gr(p), hessian = function(p) -cgradhesshalfnormAD$he(p),
        control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
          eval.max = itermax, rel.tol = tol, x.tol = tol)))
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$gradient <- cgradhesshalfnormAD$gr(mleObj$par)
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
        mleObj$hessian <- cgradhesshalfnormAD$he(mleObj$par)
      if (method == "sr1")
        mleObj$hessian <- cgradhesshalfnormAD$he(mleObj$solution)
      if (method == "mla")
        mleObj$hessian <- cgradhesshalfnormAD$he(mleObj$b)
    }
    mleObj$logL_OBS <- cgradbhhhhalfnormAD$fn(mlParam)
    mleObj$gradL_OBS <- fnGradObshalf(mlParam)
    rm(cgradbhhhhalfnormAD, cgradhesshalfnormAD, fnGradObshalf)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
    }
    ## numerical derivatives ------
    if (derivs == "numerical") {
      cgradhalfnormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::jacobian(chalfnormlike, var = unname(parm), params = list(nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar), accuracy = accuracy,
          stepsize = stepsize)
      }
      chesshalfnormlikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::hessian(function(parm) sum(chalfnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)), var = unname(parm),
          accuracy = accuracy, stepsize = stepsize)
      }
      ### solve for different algorithms ------
      mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(chalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
        gr = function(parm) -colSums(cgradhalfnormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
          maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
        maxLikAlgo = maxRoutine(fn = chalfnormlike, grad = cgradhalfnormlikeNum,
          hess = chesshalfnormlikeNum, start = startVal, finalHessian = if (hessianType ==
          2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradhalfnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), method = "SR1", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradhalfnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hs = function(parm) as(-chesshalfnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"), method = "Sparse",
          control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(chalfnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          gr = function(parm) -colSums(cgradhalfnormlikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hess = function(parm) -chesshalfnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo,
          maxiter = itermax, epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradhalfnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = function(parm) -chesshalfnormlikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
      if (method %in% c("ucminf", "nlminb")) {
        mleObj$gradient <- colSums(cgradhalfnormlikeNum(mleObj$par, nXvar = nXvar,
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
          mleObj$hessian <- chesshalfnormlikeNum(parm = mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        if (method == "sr1")
          mleObj$hessian <- chesshalfnormlikeNum(parm = mleObj$solution,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
      mleObj$logL_OBS <- chalfnormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      mleObj$gradL_OBS <- cgradhalfnormlikeNum(parm = mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
    } else {
      ## analytical derivatives ------
      if (derivs == "analytical") {
        ### solve for different algorithms ------
        mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
          fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradhalfnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
          stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = chalfnormlike,
          grad = cgradhalfnormlike, hess = chesshalfnormlike, start = startVal,
          finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradhalfnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradhalfnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hs = function(parm) as(-chesshalfnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"),
          method = "Sparse", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
          fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradhalfnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hess = function(parm) -chesshalfnormlike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo, maxiter = itermax,
          epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradhalfnormlike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = function(parm) -chesshalfnormlike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
        if (method %in% c("ucminf", "nlminb")) {
          mleObj$gradient <- colSums(cgradhalfnormlike(mleObj$par, nXvar = nXvar,
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
          mleObj$hessian <- chesshalfnormlike(parm = mleObj$par, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
          if (method == "sr1")
          mleObj$hessian <- chesshalfnormlike(parm = mleObj$solution, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        }
        mleObj$logL_OBS <- chalfnormlike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        mleObj$gradL_OBS <- cgradhalfnormlike(parm = mlParam, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
    }
  }
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
chalfnormeff <- function(object, level) {
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
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  m <- ifelse(mustar > 0, mustar, 0)
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)  ## for cost it is Farrell equivalent
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
#' marginal impact on efficiencies for halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
cmarghalfnorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}

cmarghalfnorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- crossprod(matrix(delta), t(uHvar))[1, ]
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
