################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Inefficiency structure: u_it = g(zit)u_i                                     #
#                         Modified Lee and Schmidt 1993                        #
#                          - g(zit) = exp(-eta_t * (t - T)): g(zit) = 1 for T  #
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
pgammanormlike_mols93 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, ngZGvar, gHvar,
  N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + ngZGvar)]
  eta[ngZGvar] <- 0
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -S * giepsi/gisq - exp(Wv)/(gisq * sqrt(exp(Wu)))
  sigmastar <- exp(Wv/2)/sqrt(gisq)
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mustar[i] + sigmastar[i] * qnorm(ifelse(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mustar[i]/sigmastar[i]) >=
      1, 0.9999999999, FiMat[i, ] + (1 - FiMat[i, ]) *
      pnorm(-mustar[i]/sigmastar[i]))))^(P - 1))
  }
  if (P < 0)
    return(NA)
  ll <- -1/2 * P * Wu - (TT - 1) * Wv/2 - 1/2 * log(gisq) -
    (TT - 1)/2 * log(2 * pi) - log(gamma(P)) - epsilon_isq/(2 *
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
pstgammanorm_mols93 <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, ngZGvar, gHvar, S, printInfo,
  tol, N, FiMat, wHvar, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstgammanorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGamma <- NULL
  } else {
    cat("Initialization: SFA + gamma-normal distribution...\n")
    initGamma <- maxLik::maxLik(logLik = cgammanormlike,
      start = cstgammanorm(olsObj = olsObj, epsiRes = epsiRes,
        S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
        nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]),
      grad = cgradgammanormlike, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat)
    Esti <- initGamma$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, Esti[nXvar + 3], rep(0.001, ngZGvar - 1))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "P",
    paste0("eta_", colnames(gHvar[, -ngZGvar])))
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
pgradgammanormlike_mols93 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, ngZGvar, gHvar,
  N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + ngZGvar)]
  eta[ngZGvar] <- 0
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  ewv_h <- exp(Wv/2)
  epsi_uv <- (ewv/ewu_h + S * giepsi)
  sigsta <- (ewv_h * gisq)
  musig <- epsi_uv * sqrt(gisq)/sigsta
  pmusig <- pnorm(-(musig))
  dmusig <- dnorm(-(musig), 0, 1)
  pmusig2 <- pnorm(musig)
  dmusig2 <- dnorm(musig, 0, 1)
  sigx1 <- (ewu_h * ewv_h * gisq)
  sigx2 <- (epsi_uv * ewv_h * gisq/sigsta^2)
  sigx3 <- (ewv/sigx1 - 0.5 * sigx2)
  sigx4 <- (epsi_uv/ewv_h - dmusig * sqrt(gisq)/pmusig)
  sigx5 <- (ewv * epsilon_isq/(2 * ewv)^2)
  sigx6 <- ((0.5 * (dmusig * sqrt(gisq)/pmusig) - 0.5 * (epsi_uv/ewv_h)) *
    ewv/sigx1 - 0.5 * P)
  sigx7 <- (sigx4 * sigx3 + 2 * sigx5 - 0.5 * (TT - 1))
  sigx8 <- (musig - dmusig/pmusig)
  sigZ1 <- sweep(Zigisq, MARGIN = 1, STATS = 0.5 * (epsi_uv/sqrt(gisq))/sigsta,
    FUN = "*") + sweep(Zigiepsi, MARGIN = 1, STATS = S *
    sqrt(gisq)/sigsta, FUN = "*") - sweep(Zigisq, MARGIN = 1,
    STATS = epsi_uv * ewv_h * sqrt(gisq)/sigsta^2, FUN = "*")
  sigZ2 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/gisq, FUN = "*") -
    sweep(Zigisq, MARGIN = 1, STATS = epsi_uv/gisq/gisq,
      FUN = "*")
  sigFi1 <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2,
    FUN = "*") + FiMat)
  sigFi2 <- dnorm(sigFi1, mean = 0, sd = 1)
  sigFi3 <- sweep(sigFi1, MARGIN = 1, STATS = ewv_h/sqrt(gisq),
    FUN = "*")
  sigFi4 <- sweep(sigFi3, MARGIN = 1, STATS = epsi_uv/gisq,
    FUN = "-")
  sigFi5 <- sweep((1 - FiMat)/sigFi2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  sigFi6 <- sweep((1 - FiMat)/sigFi2, MARGIN = 1, STATS = dmusig2 *
    sigx3 * sqrt(gisq), FUN = "*")
  sigFi7 <- sweep((sigFi6 + 0.5 * sigFi1), MARGIN = 1, STATS = ewv_h/sqrt(gisq),
    FUN = "*")
  sigFi8 <- sweep(sigFi7, MARGIN = 1, STATS = ewv/(ewu_h *
    gisq), FUN = "-")
  sdDiv <- apply((sigFi4)^(P - 1), 1, sum)
  sigFi9 <- sweep((sigFi5 - 1) * (sigFi4)^(P - 2), MARGIN = 1,
    STATS = S * (P - 1)/gisq, FUN = "*")
  sigFi10 <- sweep((0.5 - 0.5 * (sigFi5)) * (sigFi4)^(P - 2),
    MARGIN = 1, STATS = ewv * (P - 1)/(ewu_h * gisq), FUN = "*")
  ZZ1 <- list()
  for (i in 1:ngZGvar) {
    ZZ1[[i]] <- sweep(sigFi5, MARGIN = 1, STATS = sigZ1[,
      i] * ewv_h/sqrt(gisq), FUN = "*") - sweep(sigFi1,
      MARGIN = 1, STATS = 0.5 * (Zigisq[, i]/gisq) * ewv_h/sqrt(gisq),
      FUN = "*")
  }
  ZZ2 <- list()
  for (i in 1:ngZGvar) {
    ZZ2[[i]] <- sweep(ZZ1[[i]], MARGIN = 1, STATS = sigZ2[,
      i], FUN = "-")
  }
  ZZ3 <- list()
  for (i in 1:ngZGvar) {
    ZZ3[[i]] <- ZZ2[[i]] * (sigFi4)^(P - 2) * (P - 1)
  }
  sigFi11 <- sapply(ZZ3, FUN = function(x) apply(x, 1, sum))
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- apply(sweep(sigFi9, MARGIN = 1, STATS = Xgi[,
      k], FUN = "*"), 1, sum)/sdDiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(sigFi10, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sdDiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep((sigFi8) * (sigFi4)^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum)/sdDiv
  }
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * sigx4/sigsta,
    FUN = "*") + gx - sweep(Xepsi_i, MARGIN = 1, STATS = 1/(2 *
    ewv), FUN = "*"), gu + sweep(uHvar, MARGIN = 1, STATS = sigx6,
    FUN = "*"), gv + sweep(vHvar, MARGIN = 1, STATS = sigx7,
    FUN = "*"), apply((sigFi4)^(P - 1) * log(sigFi4), 1,
    sum)/sdDiv - (0.5 * (Wu) + digamma(P)), sweep(sigZ1,
    MARGIN = 1, STATS = sigx8, FUN = "*") + sweep(sigFi11,
    MARGIN = 1, STATS = 1/sdDiv, FUN = "*") - 0.5 * (Zigisq/gisq))
  return(sweep(gradll[, -(nXvar + nuZUvar + nvZVvar + ngZGvar +
    1)], MARGIN = 1, STATS = wHvar, FUN = "*"))
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
gammanormAlgOpt_mols93 <- function(start, olsParam, dataTable,
  S, gHvar, ngZGvar, nXvar, N, NT, FiMat_N, FiMat_NT, uHvar_c,
  uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, pindex, TT,
  Yvar, Xvar, wHvar_c, wHvar_p, method, printInfo, itermax,
  whichStart, initIter, initAlg, stepmax, tol, gradtol, hessianType,
  qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstgammanorm_mols93(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      N = NT, FiMat = FiMat_NT, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, gHvar = gHvar, ngZGvar = ngZGvar,
      uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar_c, whichStart = whichStart,
      initIter = initIter, initAlg = initAlg, tol = tol,
      printInfo = printInfo)
    initGamma <- start_st$initGamma
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(pgammanormlike_mols93(startVal, nXvar = nXvar,
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
  cat("SFA Panel MOLS93 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(pgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pgammanormlike_mols93,
    grad = pgradgammanormlike_mols93, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(pgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(pgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    hs = function(parm) as(calculus::jacobian(function(parm) -colSums(pgradgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
      parm), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(pgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    gr = function(parm) -colSums(pgradgammanormlike_mols93(parm,
      nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
    print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(pgammanormlike_mols93(parm,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)), gradient = function(parm) -colSums(pgradgammanormlike_mols93(parm,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)), control = list(iter.max = itermax,
    trace = if (printInfo) 1 else 0, eval.max = itermax,
    rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradgammanormlike_mols93(mleObj$par,
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradgammanormlike_mols93(parm,
        nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
        mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradgammanormlike_mols93(parm,
        nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)),
        mleObj$solution)
  }
  mleObj$logL_OBS <- pgammanormlike_mols93(parm = mlParam,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)
  mleObj$gradL_OBS <- pgradgammanormlike_mols93(parm = mlParam,
    nXvar = nXvar, gHvar = gHvar, ngZGvar = ngZGvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, N = N, FiMat = FiMat_N)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
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
pgammanormeff_mols93 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  eta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    2):(object$nXvar + object$nuZUvar + object$nvZVvar +
    object$ngZGvar)]
  eta[object$ngZGvar] <- 0
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
  epsilon_it <- model.response(model.frame(object$formula,
    data = object$dataTable)) - as.numeric(crossprod(matrix(beta),
    t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mui <- -(exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  Hi1 <- numeric(object$Nid)
  Hi2 <- numeric(object$Nid)
  for (i in 1:object$Nid) {
    Hi1[i] <- mean((mui[i] + sigmastar[i] * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sigmastar[i])))^(P))
    Hi2[i] <- mean((mui[i] + sigmastar[i] * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sigmastar[i])))^(P -
      1))
  }
  u <- Hi1/Hi2
  res <- data.frame(levels(pindex[, 1]), u = u, mui = mui,
    sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-u)
    mui_Gi <- res$mui - res$sigmastar^2 * git
    mui_Ki <- res$mui + res$sigmastar^2 * git
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sigmastar[i] * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sigmastar[i])))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sigmastar[i] * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sigmastar[i])))^(P -
        1))
    }
    res$teBC <- exp(-mui_Gi + sigmastar/2) * pnorm(mui_Gi/sigmastar) *
      Gi/(pnorm(mui/sigmastar) * Hi2)
    res$teBC_reciprocal <- exp(mui_Ki + sigmastar/2) * pnorm(mui_Ki/sigmastar) *
      Ki/(pnorm(mui/sigmastar) * Hi2)
  }
  res$mui <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
