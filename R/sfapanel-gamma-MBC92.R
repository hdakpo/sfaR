################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Type:      - Battese and Coelli 1992 - alternative (mbc92)                   #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = 1 + eta1 * (t - T) + eta2 * (t - T)^2                  #
# Convolution: gamma - normal                                                  #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for gamma-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
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
pgammanormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 3]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -S * giepsi/gisq - exp(Wv)/(gisq * sqrt(exp(Wu)))
  sigmastar <- exp(Wv/2)/sqrt(gisq)
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mustar[i] + sigmastar[i] * pmin(qnorm(FiMat[i, ] + (1 - FiMat[i,
      ]) * pnorm(-mustar[i]/sigmastar[i])), sqrt(.Machine$double.xmax/N)))^(P -
      1))
  }
  if (P < 0)
    return(-.Machine$double.xmax)
  ll <- -1/2 * P * Wu - (TT - 1) * Wv/2 - 1/2 * log(gisq) - (TT - 1)/2 * log(2 *
    pi) - log(gamma(P)) - epsilon_isq/(2 * exp(Wv)) + 1/2 * (mustar/sigmastar)^2 +
    pnorm(mustar/sigmastar, log.p = TRUE) + log(Hi)
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
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
pstgammanorm_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, printInfo, tol, N, wHvar, FiMat, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstgammanorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGamma <- NULL
  } else {
    cat("Initialization: SFA + gamma-normal distribution...\n")
    initGamma <- maxLik::maxLik(logLik = cgammanormlike, start = cstgammanorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradgammanormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat)
    Esti <- initGamma$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, Esti[nXvar + 3], eta1 = 0.001, eta2 = 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "P", "eta1", "eta2")
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
#' @noRd
pgradgammanormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, N, wHvar, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 3]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) rev(-(0:(x - 1)))))
  gitZit <- 2 * (Zit * (1 + Zit * (eta1 + eta2 * Zit)))
  giZi <- as.numeric(tapply(gitZit, pindex[, 1], sum))
  gitZitsq <- 2 * (Zit^2 * (1 + Zit * (eta1 + eta2 * Zit)))
  giZisq <- as.numeric(tapply(gitZitsq, pindex[, 1], sum))
  Zit_epsit <- Zit * epsilon_it
  Zi_epsi <- as.numeric(tapply(Zit_epsit, pindex[, 1], sum))
  Zitsq_epsit <- Zit^2 * epsilon_it
  Zisq_epsi <- as.numeric(tapply(Zitsq_epsit, pindex[, 1], sum))
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
  sigx6 <- ((0.5 * (dmusig * sqrt(gisq)/pmusig) - 0.5 * (epsi_uv/ewv_h)) * ewv/sigx1 -
    0.5 * P)
  sigx7 <- (sigx4 * sigx3 + 2 * sigx5 - 0.5 * (TT - 1))
  sigx8 <- (musig - dmusig/pmusig)
  sigx9 <- (0.5 * (epsi_uv * giZi/sqrt(gisq)) + S * sqrt(gisq) * Zi_epsi)
  sigx10 <- (sigx9/sigsta - epsi_uv * ewv_h * sqrt(gisq) * giZi/sigsta^2)
  sigx11 <- (S * Zi_epsi - epsi_uv * giZi/gisq)
  sigx12 <- (0.5 * (epsi_uv * giZisq/sqrt(gisq)) + S * sqrt(gisq) * Zisq_epsi)
  sigx13 <- (sigx12/sigsta - epsi_uv * ewv_h * sqrt(gisq) * giZisq/sigsta^2)
  sigx14 <- (S * Zisq_epsi - epsi_uv * giZisq/gisq)
  sigFi1 <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat)
  sigFi2 <- dnorm(sigFi1, mean = 0, sd = 1)
  sigFi3 <- sweep(sigFi1, MARGIN = 1, STATS = ewv_h/sqrt(gisq), FUN = "*")
  sigFi4 <- sweep(sigFi3, MARGIN = 1, STATS = epsi_uv/gisq, FUN = "-")
  sigFi5 <- sweep((1 - FiMat)/sigFi2, MARGIN = 1, STATS = dmusig2, FUN = "*")
  sigFi6 <- sweep((1 - FiMat)/sigFi2, MARGIN = 1, STATS = dmusig2 * sigx3 * sqrt(gisq),
    FUN = "*")
  sigFi7 <- sweep((sigFi6 + 0.5 * sigFi1), MARGIN = 1, STATS = ewv_h/sqrt(gisq),
    FUN = "*")
  sigFi8 <- sweep(sigFi7, MARGIN = 1, STATS = ewv/(ewu_h * gisq), FUN = "-")
  sdDiv <- apply((sigFi4)^(P - 1), 1, sum)
  sigFi9 <- sweep((sigFi5 - 1) * (sigFi4)^(P - 2), MARGIN = 1, STATS = S * (P -
    1)/gisq, FUN = "*")
  sigFi10 <- sweep((0.5 - 0.5 * (sigFi5)) * (sigFi4)^(P - 2), MARGIN = 1, STATS = ewv *
    (P - 1)/(ewu_h * gisq), FUN = "*")
  sigFi11 <- sweep(sigFi1, MARGIN = 1, STATS = (giZi/gisq), FUN = "*")
  sigFi12 <- sweep(sigFi5, MARGIN = 1, STATS = sigx10, FUN = "*")
  sigFi13 <- sweep((sigFi12 - 0.5 * sigFi11), MARGIN = 1, STATS = ewv_h/sqrt(gisq),
    FUN = "*")
  sigFi14 <- sweep(sigFi13, MARGIN = 1, STATS = sigx11/gisq, FUN = "-")
  sigFi15 <- sweep(sigFi1, MARGIN = 1, STATS = (giZisq/gisq), FUN = "*")
  sigFi16 <- sweep(sigFi5, MARGIN = 1, STATS = sigx13, FUN = "*")
  sigFi17 <- sweep((sigFi16 - 0.5 * sigFi15), MARGIN = 1, STATS = ewv_h/sqrt(gisq),
    FUN = "*")
  sigFi18 <- sweep(sigFi17, MARGIN = 1, STATS = sigx14/gisq, FUN = "-")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- apply(sweep(sigFi9, MARGIN = 1, STATS = Xgi[, k], FUN = "*"),
      1, sum)/sdDiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(sigFi10, MARGIN = 1, STATS = uHvar[, k], FUN = "*"),
      1, sum)/sdDiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep((sigFi8) * (sigFi4)^(P - 2) * (P - 1), MARGIN = 1,
      STATS = vHvar[, k], FUN = "*"), 1, sum)/sdDiv
  }
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * sigx4/sigsta, FUN = "*") +
    gx - sweep(Xepsi_i, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*"), gu + sweep(uHvar,
    MARGIN = 1, STATS = sigx6, FUN = "*"), gv + sweep(vHvar, MARGIN = 1, STATS = sigx7,
    FUN = "*"), apply((sigFi4)^(P - 1) * log(sigFi4), 1, sum)/sdDiv - (0.5 *
    (Wu) + digamma(P)), sigx10 * sigx8 + apply((sigFi14) * (sigFi4)^(P - 2) *
    (P - 1), 1, sum)/sdDiv - 0.5 * (giZi/gisq), sigx13 * sigx8 + apply((sigFi18) *
    (sigFi4)^(P - 2) * (P - 1), 1, sum)/sdDiv - 0.5 * (giZisq/gisq))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for gamma-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
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
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik 
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
gammanormAlgOpt_mbc92 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, N, NT, FiMat_N, FiMat_NT, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p,
  nvZVvar, Yvar, Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo, itermax,
  stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstgammanorm_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      N = NT, FiMat = FiMat_NT, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, tol = tol,
      printInfo = printInfo)
    initGamma <- start_st$initGamma
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(pgammanormlike_mbc92(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, N = N, FiMat = FiMat_N, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel Modified BC92 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(pgammanormlike_mbc92(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradgammanormlike_mbc92(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    N = N, FiMat = FiMat_N)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pgammanormlike_mbc92,
    grad = pgradgammanormlike_mbc92, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(pgammanormlike_mbc92(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), gr = function(parm) -colSums(pgradgammanormlike_mbc92(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(pgammanormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
        FiMat = FiMat_N)), gr = function(parm) -colSums(pgradgammanormlike_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(pgradgammanormlike_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), unname(parm)),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(pgammanormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N,
        FiMat = FiMat_N)), gr = function(parm) -colSums(pgradgammanormlike_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), print.info = printInfo,
      maxiter = itermax, epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(pgammanormlike_mbc92(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
        N = N, FiMat = FiMat_N)), gradient = function(parm) -colSums(pgradgammanormlike_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), control = list(iter.max = itermax,
        trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradgammanormlike_mbc92(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      N = N, FiMat = FiMat_N))
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradgammanormlike_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradgammanormlike_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- pgammanormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)
  mleObj$gradL_OBS <- pgradgammanormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, N = N, FiMat = FiMat_N)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, initGamma = initGamma))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for gamma-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
pgammanormeff_mb92 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  eta1 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 2]
  eta2 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 3]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
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
  epsilon_it <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mui <- -(exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  Hi1 <- numeric(object$Nid)
  Hi2 <- numeric(object$Nid)
  for (i in 1:object$Nid) {
    Hi1[i] <- mean((mui[i] + sigmastar[i] * qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
      ]) * pnorm(-mui[i]/sigmastar[i])))^(P))
    Hi2[i] <- mean((mui[i] + sigmastar[i] * qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
      ]) * pnorm(-mui[i]/sigmastar[i])))^(P - 1))
  }
  u <- Hi1/Hi2
  res <- data.frame(levels(pindex[, 1]), u = u, mui = mui, sigmastar = sigmastar)
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
      Gi[i] <- mean((mui_Gi[i] + sigmastar[i] * qnorm(object$FiMat[i, ] + (1 -
        object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sigmastar[i])))^(P - 1))
      Ki[i] <- mean((mui_Ki[i] + sigmastar[i] * qnorm(object$FiMat[i, ] + (1 -
        object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sigmastar[i])))^(P - 1))
    }
    res$teBC <- exp(-mui_Gi + sigmastar/2) * pnorm(mui_Gi/sigmastar) * Gi/(pnorm(mui/sigmastar) *
      Hi2)
    res$teBC_reciprocal <- exp(mui_Ki + sigmastar/2) * pnorm(mui_Ki/sigmastar) *
      Ki/(pnorm(mui/sigmastar) * Hi2)
  }
  res$mui <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
