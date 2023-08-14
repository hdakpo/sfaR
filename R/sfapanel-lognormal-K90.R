################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Type:      - Kumbakhar 1990 (K90)                                            #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = (1 + exp(eta1 * t + eta2 * t^2))^(-1)                  #
# Convolution: lognormal - normal                                              #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for lognormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
plognormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar, vHvar,
  muHvar, Yvar, Xvar, pindex, TT, S, wHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 2]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  ll <- numeric(N)
  for (i in 1:N) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(FiMat[i, ]))
    ll[i] <- -TT[i]/2 * log(2 * pi) - TT[i] * Wv[i]/2 + log(mean(exp(-1/(2 *
      exp(Wv[i])) * (epsilon_isq[i] + 2 * S * ur * giepsi[i] + gisq[i] * ur^2))))
  }
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for lognormal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
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
pstlognorm_k90 <- function(olsObj, epsiRes, nXvar, nuZUvar, muHvar, nmuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, printInfo, tol, N, FiMat, whichStart, initIter,
  initAlg) {
  if (whichStart == 1L) {
    Esti <- cstlognorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE],
      nmuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE])
    initLog <- NULL
  } else {
    cat("Initialization: SFA + lognormal-normal distribution...\n")
    initLog <- maxLik::maxLik(logLik = clognormlike, start = cstlognorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], nmuZUvar = 1, muHvar = muHvar[,
        1, drop = FALSE]), grad = cgradlognormlike, method = initAlg, control = list(iterlim = initIter,
      printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar, nuZUvar = 1,
      nmuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE], nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat)
    Esti <- initLog$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nmuZUvar > 1) {
    rep(0, nmuZUvar - 1)
  }, Esti[nXvar + 2], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 3], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_", colnames(muHvar)),
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "eta1", "eta2")
  return(list(StartVal = StartVal, initLog = initLog))
}

# Gradient of the likelihood function ----------
#' gradient for lognormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
pgradlognormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar,
  vHvar, muHvar, Yvar, Xvar, pindex, TT, S, wHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 2]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) seq(1:x)))
  gitd1 <- Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd1_epsit_gitsq <- -gitd1 * epsilon_it * git^2
  gid1_epsi_gisq <- as.numeric(tapply(gitd1_epsit_gitsq, pindex[, 1], sum))
  gitd2 <- 2 * Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd2_gitcub <- -gitd2 * git^3
  gid2_gicub <- as.numeric(tapply(gitd2_gitcub, pindex[, 1], sum))
  gitd1sq_epsit_gitsq <- -Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) * epsilon_it *
    git^2
  gid1sq_epsi_gisq <- as.numeric(tapply(gitd1sq_epsit_gitsq, pindex[, 1], sum))
  gitd2sq_gitcub <- -2 * Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) * git^3
  gid2sq_gicub <- as.numeric(tapply(gitd2sq_gitcub, pindex[, 1], sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  qFimat <- qnorm(FiMat)
  EqFi <- exp(sweep(sweep(qFimat, MARGIN = 1, STATS = ewu_h, FUN = "*"), MARGIN = 1,
    STATS = mu, FUN = "+"))
  EqFi_epsi <- sweep(sweep(EqFi, MARGIN = 1, STATS = gisq, FUN = "*"), MARGIN = 1,
    STATS = 2 * (S * giepsi), FUN = "+")
  sigFi1 <- sweep(EqFi_epsi * EqFi, MARGIN = 1, STATS = epsilon_isq, FUN = "+")
  sigFi2 <- exp(-sweep(sigFi1, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*"))
  sigDiv <- apply(sigFi2, 1, sum)
  sigFi3 <- sweep(sigFi1 * sigFi2, MARGIN = 1, STATS = ewv/(2 * ewv)^2, FUN = "*")
  sigFi4 <- sweep(sigFi2 * EqFi * qFimat, MARGIN = 1, STATS = ewu_h/(2 * ewv),
    FUN = "*") * (0.5 * EqFi_epsi + sweep(EqFi, MARGIN = 1, STATS = 0.5 * (gisq),
    FUN = "*"))
  sigFi5 <- sweep(sigFi2 * EqFi, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*") *
    (sweep(sweep(EqFi, MARGIN = 1, STATS = gid2sq_gicub, FUN = "*"), MARGIN = 1,
      STATS = 2 * (S * gid1sq_epsi_gisq), FUN = "+"))
  sigFi6 <- sweep(sigFi2 * EqFi, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*") *
    sweep(sweep(EqFi, MARGIN = 1, STATS = gid2_gicub, FUN = "*"), MARGIN = 1,
      STATS = 2 * (S * gid1_epsi_gisq), FUN = "+")
  sigFi7 <- sweep(sigFi2, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*")
  sigFi8 <- sweep(sigFi2 * EqFi, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*") *
    sweep(sweep(EqFi, MARGIN = 1, STATS = 2 * gisq, FUN = "*"), MARGIN = 1, STATS = 2 *
      (S * giepsi), FUN = "+")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- -apply(sweep(2 * (S * EqFi) * sigFi7, MARGIN = 1, STATS = Xgi[,
      k], FUN = "*") + sweep(sigFi7, MARGIN = 1, STATS = Xepsi_i[, k], FUN = "*"),
      1, sum)/sigDiv
  }
  gmu <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    gmu[, k] <- apply(sweep(sigFi8, MARGIN = 1, STATS = -muHvar[, k], FUN = "*"),
      1, sum)/sigDiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(sigFi4, MARGIN = 1, STATS = -uHvar[, k], FUN = "*"),
      1, sum)/sigDiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(2 * sigFi3, MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum)/sigDiv
  }
  gradll <- cbind(gx, gmu, gu, gv - sweep(vHvar, MARGIN = 1, STATS = 0.5 * TT,
    FUN = "*"), apply(-sigFi6, 1, sum)/sigDiv, apply(-sigFi5, 1, sum)/sigDiv)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for lognormal-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar_c matrix of Zmu variables for pooled data
#' @param uHvar_c matrix of Zu variables for pooled data
#' @param vHvar_c matrix of Zv variables for pooled data
#' @param muHvar_p matrix of Zmu variables for cross-section
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
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik 
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
lognormAlgOpt_k90 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, muHvar_c, muHvar_p, nmuZUvar, N, NT, FiMat_N, FiMat_NT, uHvar_c, uHvar_p,
  nuZUvar, vHvar_c, vHvar_p, nvZVvar, pindex, TT, Yvar, Xvar, wHvar_c, wHvar_p,
  method, printInfo, itermax, whichStart, initIter, initAlg, stepmax, tol, gradtol,
  hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstlognorm_k90(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nmuZUvar = nmuZUvar, muHvar = muHvar_c, N = NT, FiMat = FiMat_NT, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c, vHvar = vHvar_c,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol, whichStart = whichStart,
      initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    initLog <- start_st$initLog
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(plognormlike_k90(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    N = N, FiMat = FiMat_N))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel K90 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(plognormlike_k90(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradlognormlike_k90(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N, FiMat = FiMat_N)), hessian = 0,
    control = list(trace = if (printInfo) 1 else 0, maxeval = itermax, stepmax = stepmax,
      xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = plognormlike_k90,
    grad = pgradlognormlike_k90, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar,
    muHvar = muHvar_p, N = N, FiMat = FiMat_N), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(plognormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar,
      muHvar = muHvar_p, N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradlognormlike_k90(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N,
      FiMat = FiMat_N)), method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(plognormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar,
      muHvar = muHvar_p, N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradlognormlike_k90(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N,
      FiMat = FiMat_N)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(pgradlognormlike_k90(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N,
      FiMat = FiMat_N)), unname(parm)), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(plognormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar,
      muHvar = muHvar_p, N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradlognormlike_k90(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N,
      FiMat = FiMat_N)), print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(plognormlike_k90(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N, FiMat = FiMat_N)), gradient = function(parm) -colSums(pgradlognormlike_k90(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N, FiMat = FiMat_N)), control = list(iter.max = itermax,
    trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradlognormlike_k90(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N, FiMat = FiMat_N))
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradlognormlike_k90(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N,
        FiMat = FiMat_N)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradlognormlike_k90(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, N = N,
        FiMat = FiMat_N)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- plognormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    N = N, FiMat = FiMat_N)
  mleObj$gradL_OBS <- pgradlognormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    N = N, FiMat = FiMat_N)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, initLog = initLog))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# fn conditional inefficiencies ----------
#' function to estimate unconditional efficiency
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the lognormal distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon_i composite noise level
#' @param epsilon_isq composite noise square
#' @param S integer for cost/prod estimation
#' @param TT time presence vector
#' @noRd
fnCondEffLogNorm_k90 <- function(u, sigmaU, sigmaV, mu, S, TT, epsilon_isq, giepsi,
  gisq) {
  1/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) * dnorm((log(u) - mu)/sigmaU) * exp(-(epsilon_isq +
    2 * S * u * giepsi + u^2 * gisq)/(2 * sigmaV^2))
}

# fn conditional efficiencies ----------
#' function to estimate unconditional efficiency
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the lognormal distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon_i composite noise level
#' @param epsilon_isq composite noise square
#' @param S integer for cost/prod estimation
#' @param TT time presence vector
#' @noRd
fnCondBCEffLogNorm_k90 <- function(u, sigmaU, sigmaV, mu, S, TT, epsilon_isq, giepsi,
  gisq, git) {
  exp(-git * u)/u * 1/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) * dnorm((log(u) -
    mu)/sigmaU) * exp(-(epsilon_isq + 2 * S * u * giepsi + u^2 * gisq)/(2 * sigmaV^2))
}

# fn reciprocal conditional efficiencies----------
#' function to estimate unconditional efficiency
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the lognormal distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon_i composite noise level
#' @param epsilon_isq composite noise square
#' @param S integer for cost/prod estimation
#' @param TT time presence vector
#' @noRd
fnCondBCreciprocalEffLogNorm_k90 <- function(u, sigmaU, sigmaV, mu, S, TT, epsilon_isq,
  giepsi, gisq, git) {
  exp(git * u)/u * 1/(sigmaU * sigmaV^TT * (2 * pi)^(TT/2)) * dnorm((log(u) - mu)/sigmaU) *
    exp(-(epsilon_isq + 2 * S * u * giepsi + u^2 * gisq)/(2 * sigmaV^2))
}

# Conditional efficiencies estimation ----------
#' efficiencies for lognormal-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
plognormeff_k90 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  eta1 <- object$mlParam[object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1]
  eta2 <- object$mlParam[object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    2]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  pindex <- object$dataTable[, 1:2]
  invariance <- object$invariance
  if (invariance == 1) {
    muHvar_p <- apply(muHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
    uHvar_p <- apply(uHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
    vHvar_p <- apply(vHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
  } else {
    if (invariance == 2) {
      muHvar_p <- apply(muHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
    } else {
      if (invariance == 3) {
        muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
      }
    }
  }
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar_p)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar_p)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar_p)))
  TT <- as.numeric(table(pindex[, 1]))
  epsilon_it <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  u <- numeric(object$Nid)
  density_epsilon_vec <- numeric(object$Nid)
  for (i in 1:object$Nid) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i, ]))
    density_epsilon_vec[i] <- mean(1/((2 * pi)^(TT[i]/2) * exp(Wv[i]/2 * TT[i])) *
      exp(-(epsilon_isq[i] + 2 * object$S * ur * giepsi[i] + ur^2 * gisq[i])/(2 *
        exp(Wv[i]))))
    u[i] <- hcubature(f = fnCondEffLogNorm_k90, lowerLimit = 0, upperLimit = Inf,
      maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
      mu = mu[i], TT = TT[i], epsilon_isq = epsilon_isq[i], giepsi = giepsi[i],
      gisq = gisq[i], S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
  }
  res <- data.frame(levels(pindex[, 1]), u = u, giepsi = giepsi, gisq = gisq, epsilon_isq = epsilon_isq,
    density_epsilon_vec = density_epsilon_vec, Wu = Wu, Wv = Wv, mu = mu, TT = TT)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teBC <- numeric(object$Nobs)
    res$teBC_reciprocal <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      res$teBC[i] <- hcubature(f = fnCondBCEffLogNorm_k90, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(res$Wu[i]/2),
        sigmaV = exp(res$Wv[i]/2), mu = res$mu[i], TT = res$TT[i], epsilon_isq = res$epsilon_isq[i],
        giepsi = res$giepsi[i], gisq = res$gisq[i], S = object$S, git = git[i],
        vectorInterface = FALSE, tol = 1e-15)$integral/res$density_epsilon_vec[i]
      res$teBC_reciprocal[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm_k90,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(res$Wu[i]/2),
        sigmaV = exp(res$Wv[i]/2), mu = res$mu[i], TT = res$TT[i], epsilon_isq = res$epsilon_isq[i],
        giepsi = res$giepsi[i], gisq = res$gisq[i], S = object$S, git = git[i],
        vectorInterface = FALSE, tol = 1e-15)$integral/res$density_epsilon_vec[i]
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
