################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Analysis                        #
# Convolution: uniform - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf uniform-normal distribution
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
czisfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv/2)) - pnorm(S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf uniform-normal distribution
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
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
cstzisfuninorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, itermax,
  printInfo, tol) {
  cat("Initialization: SFA + uniform - normal distributions...\n")
  initUni <- maxLik(logLik = cuninormlike, start = cstuninorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgraduninormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initUni$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("SF_",
    colnames(Zvar)))
  names(initUni$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initUni = initUni))
}

# Gradient of the likelihood function ----------
#' gradient for zisf uniform-normal distribution
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
cgradzisfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  ezusq <- (sqrt(12) * (wzdeno * ewusr))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  ddue <- (depsi - dmusig)
  wzsq <- ezusq^2
  sigx1 <- (ddue * ewz/ezusq + S * prC * depsi * (epsilon)/ewvsr^2)
  sigx2 <- (prC * depsi/ewvsr + ewz * (pmusig - pepsi)/ezusq)
  sigx3 <- (dmusig/(2 * (wzdeno * ewvsr)) - sqrt(12)/2 * (wzdeno *
    ewusr * (pmusig - pepsi)/wzsq))
  sigx4 <- prC * depsi/(wzdeno * ewvsr)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * ewz/ezusq + sigx6 * prC)
  sigx8 <- (ewusr * ewz/wzsq)
  sigx9 <- (1/ezusq - sqrt(12) * sigx8)
  sigx10 <- (sigx9 * (pmusig - pepsi) - sigx4)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/(sigx2 *
    ewvsr), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx3 *
    ewz/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx7/(sigx2 *
    ewvsr), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx10 *
    ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf uniform-normal distribution
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
chesszisfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  eeusq <- (sqrt(12) * ewusr + S * (epsilon))
  ezusq <- (sqrt(12) * (wzdeno * ewusr))
  dmusig <- dnorm(eeusq/ewvsr, 0, 1)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pmusig <- pnorm(eeusq/ewvsr)
  pepsi <- pnorm(S * (epsilon)/ewvsr)
  ddue <- (depsi - dmusig)
  wzsq <- ezusq^2
  sigx1 <- (ddue * ewz/ezusq + S * prC * depsi * (epsilon)/ewvsr^2)
  sigx2 <- (prC * depsi/ewvsr + ewz * (pmusig - pepsi)/ezusq)
  sigx3 <- (dmusig/(2 * (wzdeno * ewvsr)) - sqrt(12)/2 * (wzdeno *
    ewusr * (pmusig - pepsi)/wzsq))
  sigx4 <- prC * depsi/(wzdeno * ewvsr)
  sigx5 <- (0.5 * (S * depsi * (epsilon)) - 0.5 * (eeusq *
    dmusig))
  sigx6 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 *
    depsi)
  sigx7 <- (sigx5 * ewz/ezusq + sigx6 * prC)
  sigx8 <- (ewusr * ewz/wzsq)
  sigx9 <- (1/ezusq - sqrt(12) * sigx8)
  sigx10 <- (sigx9 * (pmusig - pepsi) - sigx4)
  sigx11 <- (S^2 * (epsilon)^2/ewvsr^2 - 1)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((prC * depsi * sigx11 + ewz * (S * depsi *
      (epsilon) - eeusq * dmusig)/ezusq)/(sigx2 * ewvsr^3) -
      sigx1^2/(sigx2 * ewvsr)^2), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (eeusq * dmusig/(2 *
      (wzdeno * ewvsr^2)) - (sigx1 * sigx3/sigx2 + sqrt(12)/2 *
      (wzdeno * ddue * ewusr/wzsq))) * ewz/(sigx2 * ewvsr),
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((0.5 * (depsi * sigx11) - 0.5 * ((eeusq^2/ewvsr^2 -
    1) * dmusig)) * ewz/ezusq + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * prC * depsi * (epsilon)/ewvsr^2)/(sigx2 *
    ewvsr) - sigx7 * sigx1/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx9 * ddue - (sigx10 *
      sigx1/sigx2 + S * prC * depsi * (epsilon)/(wzdeno *
      ewvsr^2))) * ewz/(sigx2 * ewvsr), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((sqrt(12)/2 * (((sqrt(12)/2 * (dmusig/ewvsr) - 12 *
      (wzdeno^2 * ewusr * (pmusig - pepsi)/wzsq)) * ewusr +
      0.5 * (pmusig - pepsi)) * wzdeno/wzsq) + sqrt(12)/2 *
      (eeusq * dmusig)/(2 * (wzdeno * ewvsr^3))) * ewusr +
      sigx3^2 * ewz/sigx2) * ewz/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx7 * sigx3 * ewvsr/(sigx2 *
      ewvsr)^2 + (0.5 * ((sqrt(12)/2 - sqrt(12)/2 * (eeusq^2/ewvsr^2)) *
      dmusig)/(sqrt(12) * wzdeno) + sqrt(12)/2 * (sigx5 *
      wzdeno * ewusr/wzsq))/(sigx2 * ewvsr)) * ewz), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sqrt(12)/2 * (sigx9 * dmusig/ewvsr) -
      (sqrt(12)/2 * wzdeno + sqrt(12) * ((0.5 - 12 * (wzdeno^2 *
        ewusr^2/wzsq)) * ewz)) * (pmusig - pepsi)/wzsq) *
      ewusr - sigx10 * sigx3 * ewz/sigx2) * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.25 * (S^3 * depsi *
      (epsilon)^3) - 0.25 * (eeusq^3 * dmusig)) * ewz/ezusq +
      S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) -
        1) - 0.25) * prC * depsi * (epsilon)^2)/(sigx2 *
      ewvsr^3) - sigx7 * (sigx5 * ewz/ezusq + sigx6 * prC +
      0.5 * (sigx2 * ewvsr))/(sigx2 * ewvsr)^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((sigx5 * sigx9 - sigx7 * sigx10/sigx2)/ewvsr - (0.5 *
      (S^2 * depsi * (epsilon)^2/(wzdeno * ewvsr^3)) -
      0.5 * (wzdeno * depsi * ewvsr/(wzdeno * ewvsr)^2)) *
      prC) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewvsr) +
      ewvsr/(wzdeno * ewvsr)^2) * depsi - (sigx10^2/sigx2 +
      (sqrt(12) * (1 - 24 * (wzdeno * ewusr^2 * ewz/wzsq)) +
        sqrt(12)) * ewusr * (pmusig - pepsi)/wzsq)) *
      ewz + sigx9 * (pmusig - pepsi) - sigx4) * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf uniform-normal distribution
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
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
zisfuninormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfuninormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfuninormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfuninormlike, grad = cgradzisfuninormlike,
      hess = chesszisfuninormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfuninormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfuninormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfuninormlike(mleObj$par,
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
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- chesszisfuninormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfuninormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfuninormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfuninormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

# Conditional efficiencies estimation ----------

czisfuninormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u1_c1 <- -exp(Wv/2) * ((dnorm((exp(Wu) + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((exp(Wu) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) -
    object$S * epsilon
  u2_c1 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))
  u1_c2 <- rep(0, object$Nobs)
  u2_c2 <- rep(0, object$Nobs)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + theta)/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) + exp(Wv/2)))/(pnorm((theta + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv/2) + exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_c2 <- rep(1, object$Nobs)
    teBC2_c2 <- rep(1, object$Nobs)
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (pnorm((object$S * epsilon + theta)/exp(Wv/2) - exp(Wv/2)) -
        pnorm(object$S * epsilon/exp(Wv/2) - exp(Wv/2)))/(pnorm((theta +
      object$S * epsilon)/exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv/2) - exp(Wv/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv/2)))
    teBC1_reciprocal_c2 <- rep(1, object$Nobs)
    teBC2_reciprocal_c2 <- rep(1, object$Nobs)
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

czisfmarguninorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarguninorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
