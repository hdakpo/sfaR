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
# Convolution: gamma - normal                                                  #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for gamma-normal distribution
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
pgammanormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  mustar <- -S * epsilon_i/TT - exp(Wv)/(TT * sqrt(exp(Wu)))
  sigmastar <- exp(Wv/2)/sqrt(TT)
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mustar[i] + sigmastar[i] * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mustar[i]/sigmastar[i])))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  ll <- -1/2 * P * Wu - (TT - 1) * Wv/2 - 1/2 * log(TT) - (TT -
    1)/2 * log(2 * pi) - log(gamma(P)) - epsilon_isq/(2 *
    exp(Wv)) + 1/2 * (mustar/sigmastar)^2 + pnorm(mustar/sigmastar,
    log.p = TRUE) + log(Hi)
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for gamma-normal distribution
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
pstgammanorm_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, itermax, printInfo,
  tol, N, FiMat) {
  cat("Initialization: SFA + gamma-normal distribution...\n")
  initGamma <- maxLik(logLik = cgammanormlike, start = cstgammanorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradgammanormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = if (printInfo) 2 else 0, reltol = tol),
    nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
    FiMat = FiMat, )
  Esti <- initGamma$estimate
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, Esti[nXvar + 3])
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "P")
  names(initGamma$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "P")
  return(list(StartVal = StartVal, initGamma = initGamma))
}

# Gradient of the likelihood function ----------
#' gradient for gamma-normal distribution
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
pgradgammanormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
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
  ewv_h <- exp(Wv/2)
  epsi_uv <- (ewv/ewu_h + S * epsilon_i)
  musig <- epsi_uv * sqrt(TT)/(TT * ewv_h)
  pmusig <- pnorm(-(musig))
  dmusig <- dnorm(-(musig), 0, 1)
  pmusig2 <- pnorm(musig)
  dmusig2 <- dnorm(musig, 0, 1)
  sigx1 <- (TT * epsi_uv * ewv_h/(TT * ewv_h)^2)
  sigx2 <- (ewv * epsilon_isq/(2 * ewv)^2)
  ttsq <- ewv/(TT * ewu_h * ewv_h)
  sigx3 <- (ttsq - 0.5 * sigx1)
  sigx4 <- epsi_uv/TT
  sigx5 <- (epsi_uv/ewv_h - dmusig * sqrt(TT)/pmusig)
  sigFi1 <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2,
    FUN = "*") + FiMat)
  sigFi2 <- dnorm(sigFi1)
  sigFi3 <- sweep(sigFi1, MARGIN = 1, STATS = ewv_h/sqrt(TT),
    FUN = "*")
  sigFi4 <- sweep(sigFi3, MARGIN = 1, STATS = epsi_uv/TT, FUN = "-")
  sdDiv <- apply(sigFi4^(P - 1), 1, sum)
  sigFi5 <- sweep((1 - FiMat)/sigFi2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  sigFi6 <- sweep((1 - FiMat)/sigFi2, MARGIN = 1, STATS = (dmusig2 *
    sigx3) * ewv_h, FUN = "*") + sweep(0.5 * sigFi1, MARGIN = 1,
    STATS = ewv_h/sqrt(TT), FUN = "*")
  sigFi7 <- sweep(sigFi6, MARGIN = 1, STATS = (ewv/(TT * ewu_h)),
    FUN = "-")
  sigFi8 <- sweep((0.5 - 0.5 * (sigFi5)) * (sigFi4)^(P - 2),
    MARGIN = 1, STATS = ewv * (P - 1)/(TT * ewu_h), FUN = "*")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- apply(sweep(S * (sigFi5 - 1) * (sigFi4)^(P -
      2) * (P - 1), MARGIN = 1, STATS = X_iM[, k]/TT, FUN = "*"),
      1, sum)/sdDiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(sigFi8, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sdDiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(sigFi7 * (sigFi4)^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)/sdDiv
  }
  gradll <- cbind(sweep(X_iM, MARGIN = 1, STATS = S * sigx5/(TT *
    ewv_h), FUN = "*") + gx - sweep(Xepsi_i, MARGIN = 1,
    STATS = 1/ewv, FUN = "*"), gu + sweep(uHvar, MARGIN = 1,
    STATS = ((0.5 * (dmusig * sqrt(TT)/pmusig) - 0.5 * (epsi_uv/ewv_h)) *
      ttsq - 0.5 * P), FUN = "*"), gv + sweep(vHvar, MARGIN = 1,
    STATS = (sigx5 * sigx3 + 2 * sigx2 - 0.5 * (TT - 1)),
    FUN = "*"), apply((sigFi4)^(P - 1) * log(sigFi4), 1,
    sum)/sdDiv - (0.5 * (Wu) + digamma(P)))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for gamma-normal distribution
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
gammanormAlgOpt_pl81 <- function(start, olsParam, dataTable,
  S, nXvar, N, NT, FiMat_N, FiMat_NT, uHvar_c, uHvar_p, nuZUvar,
  vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c, wHvar_p,
  pindex, TT, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstgammanorm_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      N = NT, FiMat = FiMat_NT, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_c, vHvar = vHvar_c,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c,
      itermax = itermax, tol = tol, printInfo = printInfo)
    InitGamma <- start_st$initGamma
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(pgammanormlike_pl81(startVal, nXvar = nXvar,
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
    fn = function(parm) -sum(pgammanormlike_pl81(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradgammanormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pgammanormlike_pl81,
    grad = pgradgammanormlike_pl81, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N), sr1 = trust.optim(x = startVal, fn = function(parm) -sum(pgammanormlike_pl81(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N)), gr = function(parm) -colSums(pgradgammanormlike_pl81(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N)), method = "SR1", control = list(maxit = itermax,
    cgtol = gradtol, stop.trust.radius = tol, prec = tol,
    report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(pgammanormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradgammanormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), hs = function(parm) as(jacobian(function(parm) -colSums(pgradgammanormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), unname(parm)), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(pgammanormlike_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p, N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradgammanormlike_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p, N = N, FiMat = FiMat_N)), print.info = printInfo,
      maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(pgammanormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), gradient = function(parm) -colSums(pgradgammanormlike_pl81(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N)), control = list(iter.max = itermax,
      trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradgammanormlike_pl81(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(pgradgammanormlike_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p, N = N, FiMat = FiMat_N)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(pgradgammanormlike_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p, N = N, FiMat = FiMat_N)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- pgammanormlike_pl81(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)
  mleObj$gradL_OBS <- pgradgammanormlike_pl81(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
    FiMat = FiMat_N)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, InitGamma = InitGamma))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}
# Conditional efficiencies estimation ----------
#' efficiencies for gamma-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
pgammanormeff_pl81 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  u <- Hi1/Hi2
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    mui_Gi <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) -
      exp(Wv)
    mui_Ki <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) +
      exp(Wv)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv[i])))))^(P -
        1))
    }
    teBC <- exp(exp(Wv)/exp(Wu/2) + object$S * epsilon +
      exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) - object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)) * Gi/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_reciprocal <- exp(-exp(Wv)/exp(Wu/2) - object$S *
      epsilon + exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) -
      object$S * epsilon/exp(Wv/2) + exp(Wv/2)) * Ki/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    res <- bind_cols(u = u, teJLMS = teJLMS, teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- bind_cols(u = u)
  }
  return(res)
}
