################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Two types: - Pitt and Lee (1981) specification - PL81                        #
#            - Time Invariant Inefficiency                                     #
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
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
pweibullnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  if (k < 0)
    return(NA)
  ll <- numeric(N)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    ll[i] <- -TT[i]/2 * log(2 * pi) - TT[i] * Wv[i]/2 + log(mean(exp(-1/(2 *
      exp(Wv[i])) * (epsilon_isq[i] + 2 * S * ur * epsilon_i[i] +
      TT[i] * ur^2))))
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
#' @param nXvar number of main variables (inputs + env. var)
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
pstweibullnorm_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, itermax, printInfo,
  tol, N, FiMat) {
  cat("Initialization: SFA + weibull-normal distribution...\n")
  initWeibull <- maxLik(logLik = cweibullnormlike, start = cstweibullnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradweibullnormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = if (printInfo) 2 else 0, reltol = tol),
    nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, )
  Esti <- initWeibull$estimate
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, Esti[nXvar + 3])
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "k")
  names(initWeibull$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "k")
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
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
pgradweibullnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  X_iM <- apply(-Xvar, 2, function(x) tapply(x, pindex[, 1],
    sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  lFimat <- (-log(1 - FiMat))^(1/k)
  lFisq <- sweep(lFimat, MARGIN = 1, STATS = 2 * (S * ewu_h *
    epsilon_i), FUN = "*")
  lFiu <- sweep(lFimat^2, MARGIN = 1, STATS = TT * (ewu_h)^2,
    FUN = "*")
  lFiepsi <- sweep(lFisq + lFiu, MARGIN = 1, STATS = epsilon_isq,
    FUN = "+")
  lFiuv <- exp(-sweep(lFiepsi, MARGIN = 1, STATS = 1/(2 * ewv),
    FUN = "*"))
  sdFiv <- apply(lFiuv, 1, sum)
  lFi1k <- log(-log(1 - FiMat))
  lFi2k <- (-log(1 - FiMat))^(2/k)
  lFiuepsi <- sweep(lFimat, MARGIN = 1, STATS = S * epsilon_i,
    FUN = "*")
  lFi3 <- sweep(lFi2k, MARGIN = 1, STATS = TT * ewu_h, FUN = "*")
  lFi4 <- sweep(lFimat, MARGIN = 1, STATS = S * ewu_h, FUN = "*")
  lFi5 <- sweep(lFiuv, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*")
  lFi6 <- sweep((lFiuepsi + lFi3) * lFiuv, MARGIN = 1, STATS = ewu_h/(2 *
    ewv), FUN = "*")
  lFi7 <- sweep(lFi1k * (2 * (lFiuepsi) + 2 * (lFi3)) * lFiuv,
    MARGIN = 1, STATS = ewu_h/(2 * (k^2 * ewv)), FUN = "*")
  lFi8 <- sweep((lFiepsi) * lFiuv, MARGIN = 1, STATS = ewv/(2 *
    ewv)^2, FUN = "*")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- -apply(sweep(2 * lFi4 * lFi5, MARGIN = 1,
      STATS = X_iM[, k], FUN = "*") + sweep(lFi5, MARGIN = 1,
      STATS = 2 * Xepsi_i[, k], FUN = "*"), 1, sum)/sdFiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(lFi6, MARGIN = 1, STATS = -uHvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(2 * sweep(lFi8, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gradll <- cbind(gx, gu, gv - sweep(vHvar, MARGIN = 1, STATS = 0.5 *
    TT, FUN = "*"), apply(lFi7, 1, sum)/sdFiv)
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
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
weibullnormAlgOpt_pl81 <- function(start, olsParam, dataTable,
  S, nXvar, N, NT, FiMat_N, FiMat_NT, uHvar_c, uHvar_p, nuZUvar,
  vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c, wHvar_p,
  pindex, TT, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstweibullnorm_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      N = NT, FiMat = FiMat_NT, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_c, vHvar = vHvar_c,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c,
      itermax = itermax, tol = tol, printInfo = printInfo)
    InitWeibull <- start_st$initWeibull
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(pweibullnormlike_pl81(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, N = N, FiMat = FiMat_N, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(pweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pweibullnormlike_pl81,
    grad = pgradweibullnormlike_pl81, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N), sr1 = trust.optim(x = startVal, fn = function(parm) -sum(pweibullnormlike_pl81(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N)), gr = function(parm) -colSums(pgradweibullnormlike_pl81(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N)), method = "SR1", control = list(maxit = itermax,
    cgtol = gradtol, stop.trust.radius = tol, prec = tol,
    report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(pweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), hs = function(parm) as(jacobian(function(parm) -colSums(pgradweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(pweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), print.info = printInfo,
      maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(pweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), gradient = function(parm) -colSums(pgradweibullnormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), control = list(iter.max = itermax,
      trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradweibullnormlike_pl81(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N))
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
      mleObj$hessian <- jacobian(function(parm) colSums(pgradweibullnormlike_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p, N = N, FiMat = FiMat_N)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(pgradweibullnormlike_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p, N = N, FiMat = FiMat_N)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- pweibullnormlike_pl81(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N)
  mleObj$gradL_OBS <- pgradweibullnormlike_pl81(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, InitWeibull = InitWeibull))
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
fnCondEffWeibull_pl81 <- function(u, sigmaU, sigmaV, k, S, TT,
  epsilon_i, epsilon_isq) {
  u * k/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) * (u/sigmaU)^(k -
    1) * exp(-(u/sigmaU)^k) * exp(-(epsilon_isq + 2 * S *
    u * epsilon_i + u^2)/(2 * sigmaV^2))
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
fnCondBCEffWeibull_pl81 <- function(u, sigmaU, sigmaV, k, S,
  TT, epsilon_i, epsilon_isq) {
  exp(-u) * k/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) * (u/sigmaU)^(k -
    1) * exp(-(u/sigmaU)^k) * exp(-(epsilon_isq + 2 * S *
    u * epsilon_i + u^2)/(2 * sigmaV^2))
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
fnCondBCreciprocalEffWeibull_pl81 <- function(u, sigmaU, sigmaV,
  k, S, TT, epsilon_i, epsilon_isq) {
  exp(u) * k/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) * (u/sigmaU)^(k -
    1) * exp(-(u/sigmaU)^k) * exp(-(epsilon_isq + 2 * S *
    u * epsilon_i + u^2)/(2 * sigmaV^2))
}

# Conditional efficiencies estimation ----------
#' efficiencies for weibull-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
pweibullnormeff_pl81 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
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
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  u <- numeric(object$Nid)
  for (i in seq_along(1:object$Nid)) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/((2 * pi)^(TT[i]/2) * exp(Wv[i]/2 *
      TT[i])) * exp(-(epsilon_isq[i] + 2 * object$S * ur *
      epsilon_i[i] + ur^2)/(2 * exp(Wv[i]))))
    u[i] <- hcubature(f = fnCondEffWeibull_pl81, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), k = k, epsilon_i = epsilon_i[i],
      epsilon_isq = epsilon_isq[i], TT = TT[i], S = object$S,
      vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon
  }
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- numeric(object$Nid)
    teBC_reciprocal <- numeric(object$Nid)
    for (i in seq_along(1:object$Nid)) {
      ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
      density_epsilon <- mean(1/((2 * pi)^(TT[i]/2) * exp(Wv[i]/2 *
        TT[i])) * exp(-(epsilon_isq[i] + 2 * object$S *
        ur * epsilon_i[i] + ur^2)/(2 * exp(Wv[i]))))
      teBC[i] <- hcubature(f = fnCondBCEffWeibull_pl81,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        k = k, epsilon_i = epsilon_i[i], epsilon_isq = epsilon_isq[i],
        TT = TT[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
      teBC_reciprocal[i] <- hcubature(f = fnCondBCreciprocalEffWeibull_pl81,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        k = k, epsilon_i = epsilon_i[i], epsilon_isq = epsilon_isq[i],
        TT = TT[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon
    }
    res <- data.frame(levels(pindex[, 1]), u = u, teJLMS = teJLMS,
      teBC = teBC, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(levels(pindex[, 1]), u = u)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
