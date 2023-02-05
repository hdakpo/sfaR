################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Inefficiency structure: u_it = g(zit)u_i                                     #
#                         Battese and Coelli 1992 specification:               #
#                          - g(zit) = exp(-eta * (t - T))                      #
#                         Cuesta and Orea (2002), Feng and Serletis (2009)     #
#                          - g(zit) = exp(-eta1 * (t - T) - eta2 * (t - T)^2)  #
#                         Alvarez, Amsler, Orea, Schmidt (2006)                #
#                          - g(zit) = exp(eta * gHvar)                         #
#                         Kumbhakar and Wang 2005 specification:               #
#                          - g(zit) = exp(eta * (t - t1))                      #
#                         Cuesta 2000 specification:                           #
#                          - g(zit) = exp(-eta_i * (t - T))                    #
# Convolution: weibull - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for weibull-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
pweibullnormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, ngZGvar,
  gHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + ngZGvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  if (k < 0)
    return(NA)
  ll <- numeric(N)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    ll[i] <- -TT[i]/2 * log(2 * pi) - TT[i] * Wv[i]/2 + log(mean(exp(-1/(2 *
      exp(Wv[i])) * (epsilon_isq[i] + 2 * S * ur * giepsi[i] +
      gisq[i] * ur^2))))
  }
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for weibull-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param nXvar number of main variables (inputs + env. var)
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
pstweibullnorm_gzit <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, ngZGvar, gHvar, S, wHvar,
  modelType, initIter, initAlg, whichStart, printInfo, tol,
  N, FiMat) {
  if (whichStart == 1L) {
    Esti <- cstweibullnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initWeibull <- NULL
  } else {
    cat("Initialization: SFA + weibull-normal distribution...\n")
    initWeibull <- maxLik::maxLik(logLik = cweibullnormlike,
      start = cstweibullnorm(olsObj = olsObj, epsiRes = epsiRes,
        S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
        nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]),
      grad = cgradweibullnormlike, hess = chessweibullnormlike,
      method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
      )
    Esti <- initWeibull$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, Esti[nXvar + 3], if (modelType %in% c("bc92a", "kw05")) {
    0.001
  } else {
    if (modelType == "bc92b") {
      c(0.001, 0.001)
    } else {
      if (modelType == "bc92c") {
        rep(0, ngZGvar)
      } else {
        if (modelType == "c00") {
          rep(0.001, ngZGvar)
        }
      }
    }
  })
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "k",
    if (modelType %in% c("bc92a", "kw05")) {
      "eta"
    } else {
      if (modelType == "bc92b") {
        c("eta1", "eta2")
      } else {
        if (modelType == "bc92c") {
          paste0("Zg_", colnames(gHvar))
        } else {
          if (modelType == "c00") {
          paste0("eta_", colnames(gHvar))
          }
        }
      }
    })
  return(list(StartVal = StartVal, initWeibull = initWeibull))
}

# Gradient of the likelihood function ----------
#' gradient for weibull-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
pgradweibullnormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, ngZGvar,
  gHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + ngZGvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[,
    1], sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[,
    1], sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  lFimat <- (-log(1 - FiMat))^(1/k)
  lFimat2 <- (-log(1 - FiMat))^(2/k)
  lFimat3 <- log(-log(1 - FiMat))
  lFisq <- sweep(lFimat, MARGIN = 1, STATS = S * ewu_h * giepsi,
    FUN = "*")
  lFiu <- sweep(lFimat^2, MARGIN = 1, STATS = (ewu_h)^2 * gisq,
    FUN = "*")
  lFiepsi <- sweep(lFiu + 2 * (lFisq), MARGIN = 1, STATS = (epsilon_isq),
    FUN = "+")
  lFiuv <- exp(-sweep(lFiepsi, MARGIN = 1, STATS = 1/(2 * ewv),
    FUN = "*"))
  sdFiv <- apply(lFiuv, 1, sum)
  lFi1 <- sweep(lFimat, MARGIN = 1, STATS = S * ewu_h, FUN = "*")
  lFi2 <- sweep(lFiuv, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*")
  lFi3 <- sweep(lFimat2, MARGIN = 1, STATS = ewu_h * gisq,
    FUN = "*")
  lFi4 <- sweep(lFimat, MARGIN = 1, STATS = S * giepsi, FUN = "*")
  lFi5 <- sweep(((lFi3 + lFi4) * lFiuv), MARGIN = 1, STATS = ewu_h/(2 *
    ewv), FUN = "*")
  lFi6 <- sweep(lFiepsi * lFiuv, MARGIN = 1, STATS = 2 * (ewv/(2 *
    ewv)^2), FUN = "*")
  lFi7 <- sweep((2 * (lFi3) + 2 * (lFi4)) * lFiuv * lFimat3,
    MARGIN = 1, STATS = ewu_h/(2 * (k^2 * ewv)), FUN = "*")
  ZZ1 <- list()
  for (i in 1:ngZGvar) {
    ZZ1[[i]] <- sweep(lFimat^2, MARGIN = 1, STATS = (ewu_h)^2 *
      Zigisq[, i], FUN = "*")
  }
  ZZ2 <- list()
  for (i in 1:ngZGvar) {
    ZZ2[[i]] <- sweep(lFi1, MARGIN = 1, STATS = 2 * (Zigiepsi[,
      i]), FUN = "*")
  }
  ZZ3 <- list()
  for (i in 1:ngZGvar) {
    ZZ3[[i]] <- lFi2 * (ZZ1[[i]] + ZZ2[[i]])
  }
  lFi8 <- sapply(ZZ3, FUN = function(x) apply(x, 1, sum))
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- -apply(sweep(2 * (lFi1) * lFi2, MARGIN = 1,
      STATS = Xgi[, k], FUN = "*") + sweep(lFi2, MARGIN = 1,
      STATS = 2 * Xepsi_i[, k], FUN = "*"), 1, sum)/sdFiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(lFi5, MARGIN = 1, STATS = -uHvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(lFi6, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gradll <- cbind(gx, gu, gv - sweep(vHvar, MARGIN = 1, STATS = 0.5 *
    TT, FUN = "*"), apply(lFi7, 1, sum)/sdFiv, sweep(-lFi8,
    MARGIN = 1, STATS = 1/sdFiv, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for weibull-normal distribution
#' @param start starting value for optimization
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar_c matrix of Zu variables for pooled data
#' @param vHvar_c matrix of Zv variables for pooled data
#' @param uHvar_p matrix of Zu variables for cross-section
#' @param vHvar_p matrix of Zv variables for cross-section
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param modelType specification of inefficiency model G(t)u_i
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar_c vector of weights (weighted likelihood) pooled data
#' @param wHvar_p vector of weights (weighted likelihood) cross-section
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param NT number of observations (pooled data)
#' @param FiMat_NT matrix of random draws for pooled data
#' @param N number of observations (cross section)
#' @param FiMat_N matrix of random draws for cross section
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
weibullnormAlgOpt_gzit <- function(start, olsParam, dataTable,
  S, modelType, gHvar, ngZGvar, nXvar, N, NT, FiMat_N, FiMat_NT,
  uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, pindex,
  TT, Yvar, Xvar, wHvar_c, wHvar_p, method, printInfo, itermax,
  whichStart, initIter, initAlg, stepmax, tol, gradtol, hessianType,
  qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstweibullnorm_gzit(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      N = NT, FiMat = FiMat_NT, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, gHvar = gHvar, ngZGvar = ngZGvar,
      uHvar = uHvar_c, initIter = initIter, initAlg = initAlg,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S,
      modelType = modelType, wHvar = wHvar_c, whichStart = whichStart,
      tol = tol, printInfo = printInfo)
    initWeibull <- start_st$initWeibull
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(pweibullnormlike_gzit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    N = N, FiMat = FiMat_N))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel BC92-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(pweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pweibullnormlike_gzit,
    grad = pgradweibullnormlike_gzit, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(pweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(pweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    hs = function(parm) as(calculus::jacobian(function(parm) -colSums(pgradweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
      parm), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(pweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradweibullnormlike_gzit(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(pweibullnormlike_gzit(parm,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)), gradient = function(parm) -colSums(pgradweibullnormlike_gzit(parm,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)), control = list(iter.max = itermax,
    trace = if (printInfo) 1 else 0, eval.max = itermax,
    rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradweibullnormlike_gzit(mleObj$par,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N))
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradweibullnormlike_gzit(parm,
        nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
        mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradweibullnormlike_gzit(parm,
        nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
        mleObj$solution)
  }
  mleObj$logL_OBS <- pweibullnormlike_gzit(parm = mlParam,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)
  mleObj$gradL_OBS <- pgradweibullnormlike_gzit(parm = mlParam,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# fn conditional inefficiencies ----------
#' function to estimate unconditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon_i composite noise level
#' @param epsilon_isq composite noise square
#' @param S integer for cost/prod estimation
#' @param TT time presence vector
#' @noRd
fnCondEffWeibull_gzit <- function(u, sigmaU, sigmaV, k, S, TT,
  epsilon_isq, giepsi, gisq) {
  u * k/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) * (u/sigmaU)^(k -
    1) * exp(-(u/sigmaU)^k) * exp(-(epsilon_isq + 2 * S *
    u * giepsi + u^2 * gisq)/(2 * sigmaV^2))
}

# fn conditional efficiencies ----------
#' function to estimate unconditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon_i composite noise level
#' @param epsilon_isq composite noise square
#' @param S integer for cost/prod estimation
#' @param TT time presence vector
#' @noRd
fnCondBCEffWeibull_gzit <- function(u, sigmaU, sigmaV, k, S,
  TT, epsilon_isq, giepsi, gisq, git) {
  exp(-git * u) * k/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) *
    (u/sigmaU)^(k - 1) * exp(-(u/sigmaU)^k) * exp(-(epsilon_isq +
    2 * S * u * giepsi + u^2 * gisq)/(2 * sigmaV^2))
}

# fn reciprocal conditional efficiencies----------
#' function to estimate unconditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon_i composite noise level
#' @param epsilon_isq composite noise square
#' @param S integer for cost/prod estimation
#' @param TT time presence vector
#' @noRd
fnCondBCreciprocalEffWeibull_gzit <- function(u, sigmaU, sigmaV,
  k, S, TT, epsilon_isq, giepsi, gisq, git) {
  exp(git * u) * k/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) *
    (u/sigmaU)^(k - 1) * exp(-(u/sigmaU)^k) * exp(-(epsilon_isq +
    2 * S * u * giepsi + u^2 * gisq)/(2 * sigmaV^2))
}

# Conditional efficiencies estimation ----------
#' efficiencies for weibull-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
pweibullnormeff_gzit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  eta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + object$nvZVvar +
    object$ngZGvar + 1)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  gHvar <- object$gHvar
  pindex <- object$dataTable[, 1:2]
  invariance <- object$invariance
  if (invariance == 1) {
    uHvar_p <- apply(uHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
    vHvar_p <- apply(vHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
  } else {
    if (invariance == 2) {
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
    } else {
      if (invariance == 3) {
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
      }
    }
  }
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar_p)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar_p)))
  TT <- as.numeric(table(pindex[, 1]))
  epsilon_it <- model.response(model.frame(object$formula,
    data = object$dataTable)) - as.numeric(crossprod(matrix(beta),
    t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  u <- numeric(object$Nid)
  density_epsilon_vec <- numeric(object$Nid)
  for (i in seq_along(1:object$Nid)) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon_vec[i] <- mean(1/((2 * pi)^(TT[i]/2) *
      exp(Wv[i]/2 * TT[i])) * exp(-(epsilon_isq[i] + 2 *
      object$S * ur * giepsi[i] + ur^2 * gisq[i])/(2 *
      exp(Wv[i]))))
    u[i] <- hcubature(f = fnCondEffWeibull_gzit, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), k = k, epsilon_isq = epsilon_isq[i],
      giepsi = giepsi[i], gisq = gisq[i], TT = TT[i], S = object$S,
      vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
  }
  res <- data.frame(levels(pindex[, 1]), u = u, giepsi = giepsi,
    gisq = gisq, epsilon_isq = epsilon_isq, density_epsilon_vec = density_epsilon_vec,
    Wu = Wu, Wv = Wv, mu = mu, TT = TT)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teBC <- numeric(object$Nid)
    res$teBC_reciprocal <- numeric(object$Nid)
    for (i in seq_along(1:object$Nid)) {
      res$teBC[i] <- hcubature(f = fnCondBCEffWeibull_gzit,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(res$Wu[i]/2), sigmaV = exp(res$Wv[i]/2),
        k = k, epsilon_isq = res$epsilon_isq[i], giepsi = res$giepsi[i],
        gisq = res$gisq[i], TT = res$TT[i], S = object$S,
        git = git[i], vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
      res$teBC_reciprocal[i] <- hcubature(f = fnCondBCreciprocalEffWeibull_gzit,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(res$Wu[i]/2), sigmaV = exp(res$Wv[i]/2),
        k = k, epsilon_isq = res$epsilon_isq[i], giepsi = res$giepsi[i],
        gisq = res$gisq[i], TT = res$TT[i], S = object$S,
        git = git[i], vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
    }
  }
  res$giepsi <- NULL
  res$gisq <- NULL
  res$epsilon_isq <- NULL
  res$density_epsilon_vec <- NULL
  res$Wu <- NULL
  res$Wv <- NULL
  res$mu <- NULL
  res$TT <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
