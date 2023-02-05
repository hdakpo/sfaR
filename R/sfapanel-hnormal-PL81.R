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
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phalfnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  mustar <- -exp(Wu) * S * epsilon_i/(exp(Wv) + TT * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + TT * exp(Wu)))
  ll <- log(2) - TT/2 * log(2 * pi) - (TT - 1)/2 * Wv - 1/2 *
    log(exp(Wv) + TT * exp(Wu)) + pnorm(mustar/sigmastar,
    log.p = TRUE) - epsilon_isq/(2 * exp(Wv)) + 1/2 * (mustar/sigmastar)^2
  return(wHvar * ll)
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
psthalfnorm_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, whichStart,
  initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA + halfnormal - normal distributions...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[,
        1]), vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  })
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
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
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
pgradhalfnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  X_iM <- apply(-Xvar, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  sigmasq <- (ewv + TT * ewu)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  musig <- (S * ewu * epsilon_i/(ssqx2))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig)
  dpmu <- dmusig/pmusig
  sigx1 <- (S * ewu * epsilon_i/(ssqx2) - dpmu)
  wvsq <- (1 - ewv/sigmasq)
  sigx2 <- (wvsq * ewu/sigmastar)
  wu2sq <- (1 - TT * ewu/sigmasq)
  sigx3 <- (dpmu - S * ewu * epsilon_i/(ssqx2))
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + TT * sigmastar)
  sigx5 <- (1/(ssqx2) - sigx4 * ewu/ssq)
  pmustar <- pmusig * sigmastar
  ttsq <- (TT/sigmasq)
  sigx6 <- (0.5 * sigx2 + sigmastar)
  epsivsq <- (epsilon_isq/(2 * ewv)^2)
  epsiusq <- ewu * epsilon_i/ssq
  gradll <- cbind(sweep(X_iM, MARGIN = 1, STATS = S * ewu *
    (S * epsilon_i/ewv - dmusig/(pmustar))/sigmasq, FUN = "*") -
    sweep(Xepsi_i, MARGIN = 1, STATS = 2/(2 * ewv), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ewu * (S * sigx5 * sigx1 *
      epsilon_i - 0.5 * ttsq), FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = ((2 * epsivsq + S * sigx6 * sigx3 *
        epsiusq - 0.5/sigmasq) * ewv - 0.5 * (TT - 1)),
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
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
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phesshalfnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  X_iM <- apply(-Xvar, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) {
      tapply(x, pindex[, 1], sum)
    })
  }
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  sigmasq <- (ewv + TT * ewu)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  musig <- (S * ewu * epsilon_i/(ssqx2))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig)
  dpmu <- dmusig/pmusig
  sigx1 <- (S * ewu * epsilon_i/(ssqx2) - dpmu)
  wvsq <- (1 - ewv/sigmasq)
  sigx2 <- (wvsq * ewu/sigmastar)
  wu2sq <- (1 - TT * ewu/sigmasq)
  sigx3 <- (dpmu - S * ewu * epsilon_i/(ssqx2))
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + TT * sigmastar)
  sigx5 <- (1/(ssqx2) - sigx4 * ewu/ssq)
  pmustar <- pmusig * sigmastar
  ttsq <- (TT/sigmasq)
  sigx6 <- (0.5 * sigx2 + sigmastar)
  epsivsq <- (epsilon_isq/(2 * ewv)^2)
  ttusq <- TT * ewu/sigmasq
  sigx7 <- dmusig * (dmusig/(pmustar) - S * epsilon_i/ewv)/pmusig
  wusq <- ewu/sigmasq
  epsiusq <- ewu * epsilon_i/ssq
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar +
    nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(X_iM, MARGIN = 1,
    STATS = S^2 * wHvar * (1/ewv - dmusig * (dmusig/(pmustar)^2 -
      S * epsilon_i/(ewv * pmustar)) * wusq) * wusq, FUN = "*"),
    X_iM) - sapply(1:nXvar, function(x) {
    crossprod(Xsq[[x]], as.matrix(wHvar * 2/(2 * ewv)))
  })
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(X_iM,
    MARGIN = 1, STATS = S * sigx5 * ewu * wHvar * (S * (2/sigmastar -
      sigx7) * ewu * epsilon_i/sigmasq - dpmu), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xepsi_i, MARGIN = 1, STATS = 2 *
    wHvar * (2/(2 * ewv)^2) * ewv, FUN = "*") + sweep(X_iM,
    MARGIN = 1, STATS = wHvar * (S * sigx6 * (dpmu + S *
      (sigx7 - 2/sigmastar) * ewu * epsilon_i/sigmasq) *
      ewu/ssq) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (TT^2/sigmasq^2) + S * (S * (1 - dmusig * sigx3/pmusig) *
      sigx5^2 * epsilon_i - ((0.5 * (ttusq) + 1 - 0.5 *
      (0.5 * wu2sq + ttusq)) * wu2sq * ewv/sigmastar +
      (2 * TT - 2 * (sigx4^2 * ewu * sigmasq/ssq)) * sigmastar) *
      sigx1/ssq) * epsilon_i) * ewu + S * sigx5 * sigx1 *
      epsilon_i - 0.5 * ttsq) * ewu, FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.5/sigmasq^2 - 16 * (ewv *
      epsilon_isq/(2 * ewv)^4)) * ewv + 2 * epsivsq + S *
      (((0.5 * (ewv/sigmasq) - 0.5 * (0.5 * wvsq + ewv/sigmasq)) *
        wvsq * sigx3/sigmastar + S * sigx6^2 * (1 + dmusig *
        sigx1/pmusig) * ewv * epsilon_i/ssq) * ewu +
        sigx6 * (1 - 2 * (sigx6 * ewv * ssqx2/ssq)) *
          sigx3) * epsiusq - 0.5/sigmasq) * ewv, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (0.5 * (TT/sigmasq^2) + S *
      (((0.5 * (wu2sq * ewv/sigmasq) + 0.5 * (1 + ewv *
        (ttusq - 1)/sigmasq - 0.5 * (wvsq * wu2sq))) *
        sigx3/sigmastar + S * sigx6 * sigx5 * (dmusig *
        sigx3/pmusig - 1) * epsilon_i) * ewu + sigx6 *
        (1 - 2 * (sigx4 * ewu * ssqx2/ssq)) * sigx3) *
      epsilon_i/ssq) * ewu * ewv, FUN = "*"), vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for halfnormal-normal distribution
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
halfnormAlgOpt_pl81 <- function(start, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar,
  Yvar, Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo,
  whichStart, initIter, initAlg, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psthalfnorm_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar_c, tol = tol, whichStart = whichStart,
      initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(phalfnormlike_pl81(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel PL81 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = phalfnormlike_pl81,
    grad = pgradhalfnormlike_pl81, hess = phesshalfnormlike_pl81,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesshalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesshalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(phalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradhalfnormlike_pl81(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phesshalfnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradhalfnormlike_pl81(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
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
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$hessian <- phesshalfnormlike_pl81(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesshalfnormlike_pl81(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- phalfnormlike_pl81(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradhalfnormlike_pl81(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
phalfnormeff_pl81 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon_i/(exp(Wv) + TT *
    exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + TT * exp(Wu)))
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
    res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB,
      uUB = uUB, teJLMS = teJLMS, m = m, teMO = teMO, teBC = teBC,
      teBCLB = teBCLB, teBCUB = teBCUB, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB,
      uUB = uUB, m = m)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
