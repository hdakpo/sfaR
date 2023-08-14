################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Type:      - Pitt and Lee (1981) specification - PL81                        #
#            - Time Invariant Inefficiency                                     #
# Convolution: rayleigh - normal                                               #
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
praynormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  mustar <- -exp(Wu) * S * epsilon_i/(exp(Wv) + TT * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + TT * exp(Wu)))
  ll <- 1/2 * (mustar/sigmastar)^2 - epsilon_isq/(2 * exp(Wv)) - (TT - 1)/2 * log(2 *
    pi) - (TT - 1)/2 * Wv - 1/2 * Wu - 1/2 * log(exp(Wv) + TT * exp(Wu)) + log(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
pstraynorm_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- cstraynorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initRay <- NULL
  } else {
    cat("Initialization: SFA + rayleigh-normal distribution...\n")
    initRay <- maxLik::maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradraynormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initRay$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  })
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)))
  return(list(StartVal = StartVal, initRay = initRay))
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
pgradraynormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
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
  wvsq <- (1 - ewv/sigmasq)
  wu2sq <- (1 - TT * ewu/sigmasq)
  squv <- (sigmastar/ewv - ewu/(ssqx2))
  wvdsig <- (wvsq * dmusig/sigmastar)
  dmuvsig <- (dmusig * ewv/sigmastar)
  sigx1 <- (dmusig * sigmastar - S * ewu * pmusig * epsilon_i/sigmasq)
  sigx2 <- (pmusig + S * dmusig * squv * epsilon_i)
  sigx3 <- (S * epsilon_i/ewv - sigx2/sigx1)
  sigx4 <- (0.5 * dmuvsig - S * pmusig * epsilon_i)
  sigx5 <- (0.5 * (wu2sq * ewv/sigmastar) + TT * sigmastar)
  sigx6 <- (1/(ssqx2) - sigx5 * ewu/ssq)
  sigx7 <- (sigx4 * wu2sq/sigx1 + S^2 * sigx6 * ewu * epsilon_i^2/sigmastar - 0.5 *
    TT)
  sigx8 <- (0.5 * wvdsig + S * pmusig * epsilon_i/sigmasq)
  sigx9 <- (0.5 * (wvsq * ewu/sigmastar) + sigmastar)
  sigx10 <- (sigx8/sigx1 - S^2 * sigx9 * ewu * epsilon_i^2/(ssq * sigmastar))
  sigx11 <- (sigx10 * ewu - 0.5)/sigmasq
  gradll <- cbind(sweep(X_iM, MARGIN = 1, STATS = S * ewu * sigx3/sigmasq, FUN = "*") -
    sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = (sigx7 * ewu/sigmasq - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = ((sigx11 + 2 * (epsilon_isq/(2 * ewv)^2)) * ewv - 0.5 * (TT - 1)),
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
phessraynormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  X_iM <- apply(-Xvar, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) {
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
  wvsq <- (1 - ewv/sigmasq)
  wu2sq <- (1 - TT * ewu/sigmasq)
  squv <- (sigmastar/ewv - ewu/(ssqx2))
  wvdsig <- (wvsq * dmusig/sigmastar)
  dmuvsig <- (dmusig * ewv/sigmastar)
  sigx1 <- (dmusig * sigmastar - S * ewu * pmusig * epsilon_i/sigmasq)
  sigx2 <- (pmusig + S * dmusig * squv * epsilon_i)
  sigx3 <- (S * epsilon_i/ewv - sigx2/sigx1)
  sigx4 <- (0.5 * dmuvsig - S * pmusig * epsilon_i)
  sigx5 <- (0.5 * (wu2sq * ewv/sigmastar) + TT * sigmastar)
  sigx6 <- (1/(ssqx2) - sigx5 * ewu/ssq)
  sigx7 <- (sigx4 * wu2sq/sigx1 + S^2 * sigx6 * ewu * epsilon_i^2/sigmastar - 0.5 *
    TT)
  sigx8 <- (0.5 * wvdsig + S * pmusig * epsilon_i/sigmasq)
  sigx9 <- (0.5 * (wvsq * ewu/sigmastar) + sigmastar)
  sigx10 <- (sigx8/sigx1 - S^2 * sigx9 * ewu * epsilon_i^2/(ssq * sigmastar))
  sigx11 <- (sigx10 * ewu - 0.5)/sigmasq
  sigx12 <- dmusig * epsilon_i/sigmastar
  sigx13 <- (TT * ewu/sigmasq)
  sigx14 <- (ssq * sigmastar)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar + nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(X_iM, MARGIN = 1, STATS = S^2 * wHvar *
    (1/ewv - (((1 - S^2 * ewu * epsilon_i^2/(ewv * sigmasq)) * squv - ewu/(ssqx2)) *
      dmusig + ewu * sigx2^2/(sigx1 * sigmasq))/sigx1) * ewu/sigmasq, FUN = "*"),
    X_iM) - sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar/ewv)))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(X_iM, MARGIN = 1,
    STATS = S * wHvar * (((sigx4 * sigx2/sigx1 + 0.5 * (S * sigx12)) * ewu/sigmasq -
      pmusig) * wu2sq/sigx1 + 2 * (S * sigx6 * ewu * epsilon_i/sigmastar)) *
      ewu/sigmasq, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xepsi_i,
    MARGIN = 1, STATS = wHvar * 2 * (2/(2 * ewv)^2) * ewv, FUN = "*") + sweep(X_iM,
    MARGIN = 1, STATS = S * wHvar * (((sigx8 * sigx2/sigx1 - S * (0.5 * (wvsq/ewv) +
      1/sigmasq) * sigx12) * ewu + pmusig)/(sigx1 * sigmasq) - 2 * (S * sigx9 *
      ewu * epsilon_i/sigx14)) * ewu/sigmasq * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * dmuvsig + ewu * (S^2 * sigx6 * dmusig *
      epsilon_i^2 - (sigx4 * wu2sq/sigx1 + TT) * sigx4/sigmasq) - (0.5 * (0.5 *
      (wu2sq * dmusig * ewv/sigmastar) + S^2 * sigx6 * dmusig * ewu * epsilon_i^2) +
      S * pmusig * epsilon_i)) * wu2sq/sigx1 + ewu * (S^2 * (2/(ssqx2) - (((0.5 *
      sigx13 + 2 - 0.5 * (0.5 * wu2sq + TT * ewu/sigmasq)) * wu2sq * ewv/sigmastar +
      (4 * TT - 2 * (sigx5^2 * ewu * sigmasq/ssq)) * sigmastar) * ewu/ssq +
      0.5 * (wu2sq * sigx6))) * epsilon_i^2/sigmastar - TT * sigx7/sigmasq) -
      0.5 * TT) * ewu/sigmasq, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((0.5 *
    (wvsq * dmusig) + 0.5 * (ewu * (TT * dmusig * ewv/sigmasq - S^2 * wvsq *
    sigx6 * dmusig * ewu * epsilon_i^2/sigmastar)/sigmasq - 0.5 * (wvsq * wu2sq *
    dmusig)))/sigmastar + (S * pmusig * epsilon_i - (sigx8 * sigx4 * wu2sq/sigx1 +
    S * (S * sigx6 * dmusig * epsilon_i + TT * pmusig/sigmasq) * epsilon_i) *
    ewu)/sigmasq)/sigx1 - (S^2 * (((0.5 * (wu2sq * ewv/sigmasq) + 0.5 * wvsq +
    0.5 * (1 + ewv * (TT * ewu/sigmasq - 1)/sigmasq - 0.5 * (wvsq * wu2sq))) *
    ewu/sigmastar + sigmastar)/sigx14 + sigx9 * (1/sigx14 - (0.5 * (ssq * wu2sq/(ssqx2)) +
    2 * (sigx5 * ewu)) * ewu * ewv/sigx14^2)) * ewu * epsilon_i^2 + TT * sigx11)) *
    ewu * ewv/sigmasq, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((0.5 * (ewv * (S^2 * sigx9 * dmusig * ewu^2 * epsilon_i^2/sigx14 - dmusig)/sigmasq -
      0.5 * (wvsq * dmusig)) + 0.5 * dmusig) * wvsq/sigmastar + (ewv * (S *
      (S * sigx9 * dmusig * ewu * epsilon_i/ssq - pmusig/sigmasq) * epsilon_i -
      sigx8^2 * ewu/sigx1) + S * pmusig * epsilon_i)/sigmasq)/sigx1 - S^2 *
      (sigx9 * (1/sigx14 - (0.5 * (ssq * wvsq/(ssqx2)) + 2 * (sigx9 * ewv)) *
        ewu * ewv/sigx14^2) + (0.5 * (ewv/sigmasq) - 0.5 * (0.5 * wvsq +
        ewv/sigmasq)) * wvsq * sigmasq/(ssq * ewv)) * ewu * epsilon_i^2) *
      ewu - ((sigx10 * ewu - 0.5) * ewv/sigmasq + 0.5))/sigmasq + (2 - 16 *
      (ewv^2/(2 * ewv)^2)) * epsilon_isq/(2 * ewv)^2) * ewv, FUN = "*"), vHvar)
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
raynormAlgOpt_pl81 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, itermax, stepmax, tol, gradtol, whichStart,
  initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstraynorm_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    InitRay <- start_st$initRay
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(praynormlike_pl81(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel PL81 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(praynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = praynormlike_pl81,
    grad = pgradraynormlike_pl81, hess = phessraynormlike_pl81, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(praynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(praynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(praynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(praynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessraynormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradraynormlike_pl81(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
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
      mleObj$hessian <- phessraynormlike_pl81(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessraynormlike_pl81(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- praynormlike_pl81(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradraynormlike_pl81(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, InitRay = InitRay))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for rayleigh-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
praynormeff_pl81 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
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
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  mustar <- -exp(Wu) * object$S * epsilon_i/(exp(Wv) + TT * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + TT * exp(Wu)))
  u <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  m <- ifelse(mustar/2 + sqrt(sigmastar^2 + mustar^2/4) > 0, mustar/2 + sqrt(sigmastar^2 +
    mustar^2/4), 0)
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar - sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teBC_reciprocal <- exp(mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar +
      sigmastar) + (mustar + sigmastar^2) * pnorm(mustar/sigmastar + sigmastar))/(sigmastar *
      dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
    teMO <- exp(-m)
    res <- data.frame(levels(pindex[, 1]), u = u, teJLMS = teJLMS, teBC = teBC,
      m = m, teMO = teMO, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(levels(pindex[, 1]), u = u, m = m)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
