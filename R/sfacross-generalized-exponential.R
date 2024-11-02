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
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  ll <- (log(2) - 1/2 * Wu + log(exp(A) * pnorm(a) - exp(B) *
    pnorm(b)))
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
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  .e1 <- Wv
  .e2 <- Wu
  .e3 <- exp(.e1/2)
  .e5 <- exp(.e2/2)
  .e7 <- S * epsilon
  .e8 <- .e3/.e5
  .e9 <- .e7/.e3
  .e10 <- exp(.e2)
  .e11 <- exp(.e1)
  .e12 <- .e7/.e5
  .e13 <- -(2 * .e8 + .e9)
  .e14 <- -(.e8 + .e9)
  .e15 <- 2 * .e10
  .e16 <- 2 * (.e11/.e10)
  .e18 <- exp(.e16 + 2 * .e12)
  .e19 <- exp(.e11/.e15 + .e12)
  .e20 <- pnorm(.e13)
  .e21 <- pnorm(.e14)
  .e22 <- dnorm(.e13, 0, 1)
  .e23 <- dnorm(.e14, 0, 1)
  .e26 <- .e19 * .e21 - .e18 * .e20
  .e27 <- 0.5 * .e9
  gradll <- cbind(S * Xvar * (((.e23/.e3 - .e21/.e5) * .e19 - (.e22/.e3 - 2 * (.e20/.e5)) *
    .e18)/.e26), uHvar * (((0.5 * (.e23 * .e3/.e5) - (0.5 * .e12 + 2 * (.e10 *
    .e11/.e15^2)) * .e21) * .e19 - (.e22 * .e3/.e5 - (.e16 + .e12) * .e20) *
    .e18)/.e26 - 0.5), vHvar * ((.e19 * (.e11 * .e21/.e15 - (0.5 * .e8 - .e27) *
    .e23) - (2 * (.e11 * .e20/.e10) - .e22 * (.e8 - .e27)) * .e18)/.e26))
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
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  .e1 <- Wv
  .e2 <- Wu
  .e3 <- exp(.e1/2)
  .e5 <- exp(.e2/2)
  .e7 <- S * epsilon
  .e8 <- .e3/.e5
  .e9 <- .e7/.e3
  .e10 <- exp(.e2)
  .e11 <- exp(.e1)
  .e12 <- .e7/.e5
  .e14 <- 2 * .e8 + .e9
  .e15 <- .e8 + .e9
  .e16 <- -.e14
  .e17 <- -.e15
  .e18 <- 2 * .e10
  .e19 <- 2 * (.e11/.e10)
  .e20 <- pnorm(.e16)
  .e21 <- pnorm(.e17)
  .e22 <- .e11/.e18
  .e24 <- exp(.e19 + 2 * .e12)
  .e25 <- exp(.e22 + .e12)
  .e26 <- dnorm(.e16, 0, 1)
  .e27 <- dnorm(.e17, 0, 1)
  .e28 <- 0.5 * .e9
  .e29 <- .e18^2
  .e30 <- 0.5 * .e8
  .e33 <- .e25 * .e21 - .e24 * .e20
  .e34 <- .e30 - .e28
  .e35 <- .e8 - .e28
  .e36 <- 0.5 * .e12
  .e37 <- 2 * (.e10 * .e11/.e29)
  .e38 <- .e36 + .e37
  .e39 <- .e19 + .e12
  .e41 <- .e26 * .e3/.e5
  .e42 <- .e34 * .e27
  .e43 <- .e26 * .e35
  .e44 <- .e27 * .e3
  .e46 <- 0.5 * (.e44/.e5)
  .e48 <- 2 * (.e11 * .e20/.e10) - .e43
  .e50 <- .e41 - .e39 * .e20
  .e52 <- .e26/.e3 - 2 * (.e20/.e5)
  .e54 <- .e27/.e3 - .e21/.e5
  .e57 <- .e11 * .e21/.e18 - .e42
  .e59 <- .e46 - .e38 * .e21
  .e61 <- .e59 * .e25 - .e50 * .e24
  .e65 <- .e54 * .e25 - .e52 * .e24
  .e67 <- .e25 * .e57 - .e48 * .e24
  .e68 <- .e34 * .e15
  .e69 <- .e14 * .e35
  .e72 <- .e61 * .e65/.e33
  .e74 <- .e61 * .e67/.e33
  .e76 <- .e65 * .e67/.e33
  .e77 <- .e68 + 0.5
  .e78 <- .e69 + 0.5
  .e79 <- 0.25 * .e9
  .e80 <- 2 * .e50
  .e81 <- 2 * (.e10 * .e21/.e29)
  hessll <- matrix(0, nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar + nuZUvar +
    nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * ((((.e15/.e3 - 1/.e5) * .e27/.e3 -
    .e54/.e5) * .e25 - (((.e14/.e3 - 2/.e5) * .e26/.e3 - 2 * (.e52/.e5)) * .e24 +
    .e65^2/.e33))/.e33) * wHvar, Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(S * Xvar * (((((0.5 +
    .e36 + .e37) * .e21 - .e46)/.e5 + (0.5 * (.e15/.e5) - .e38/.e3) * .e27) *
    .e25 - (((.e14/.e5 - .e39/.e3) * .e26 + (.e20 - .e80)/.e5) * .e24 + .e72))/.e33) *
    wHvar, uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(S *
    Xvar * (((.e27 * (.e22 - .e77)/.e3 - .e57/.e5) * .e25 - (((.e19 - .e78) *
    .e26/.e3 - 2 * (.e48/.e5)) * .e24 + .e76))/.e33) * wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(uHvar *
    ((((0.5 * (0.5 * (.e3 * .e15/.e5) - 0.5) - 0.5 * .e38) * .e27 * .e3/.e5 -
      (.e59 * .e38 + (2 * ((1 - 8 * (.e10^2/.e29)) * .e10 * .e11/.e29) - 0.25 *
        .e12) * .e21)) * .e25 - ((((.e14 * .e3 - .e7)/.e5 - (0.5 + .e19)) *
      .e26 * .e3/.e5 + (.e36 + .e19) * .e20 - .e39 * .e50) * .e24 + .e61^2/.e33))/.e33) *
    wHvar, uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(uHvar * ((((.e44/(4 * (.e10 * .e5)) - .e81) * .e11 -
    ((0.5 * .e68 - 0.25) * .e27 * .e3/.e5 + .e38 * .e57)) * .e25 - (.e74 + (2 *
    ((.e41 - .e20) * .e11/.e10) - ((.e69 - 0.5) * .e26 * .e3/.e5 + .e48 * .e39)) *
    .e24))/.e33) * wHvar, vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(vHvar * ((((.e57/2 + (.e21 -
    .e42)/2) * .e11/.e10 - (0.25 * .e8 + .e79 - .e34^2 * .e15) * .e27) * .e25 -
    (((2 * .e48 + 2 * (.e20 - .e43)) * .e11/.e10 - (.e79 + .e30 - .e14 * .e35^2) *
      .e26) * .e24 + .e67^2/.e33))/.e33) * wHvar, vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for generalized exponential-normal distribution
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
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
genexponormAlgOpt <- function(start, olsParam, dataTable, S,
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar,
  method, printInfo, itermax, stepmax, tol, gradtol, hessianType,
  qac) {
  startVal <- if (!is.null(start))
    start else cstgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(cgenexponormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cgenexponormlike,
    grad = cgradgenexponormlike, hess = chessgenexponormlike,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S)), hs = function(parm) as(-chessgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S)), hess = function(parm) -chessgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S)), hessian = function(parm) -chessgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S), control = list(iter.max = itermax,
        trace = if (printInfo) 1 else 0, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradgenexponormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S))
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
      mleObj$hessian <- chessgenexponormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S)
    if (method == "sr1")
      mleObj$hessian <- chessgenexponormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S)
  }
  mleObj$logL_OBS <- cgenexponormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S)
  mleObj$gradL_OBS <- cgradgenexponormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
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
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
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
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
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
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
