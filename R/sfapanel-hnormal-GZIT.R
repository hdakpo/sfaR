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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phalfnormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, ngZGvar, gHvar,
  wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + ngZGvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -exp(Wu) * S * giepsi/(exp(Wv) + gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  ll <- log(2) - TT/2 * log(2 * pi) + log(pnorm(mustar/sigmastar)) +
    log(sigmastar) - 1/2 * (epsilon_isq/exp(Wv) - (mustar/sigmastar)^2) -
    TT/2 * log(exp(Wv)) - 1/2 * log(exp(Wu))
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
#' @param modelType specification of inefficiency model G(t)u_i
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param whichStart strategy to get starting values
#' @param printInfo logical print info during optimization
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
psthalfnorm_gzit <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, ngZGvar, gHvar, Yvar, Xvar, S, wHvar,
  modelType, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA + halfnormal-normal distribution...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, if (modelType %in% c("bc92a", "kw05")) {
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
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), if (modelType %in%
    c("bc92a", "kw05")) {
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
pgradhalfnormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, ngZGvar, gHvar,
  wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + ngZGvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[,
    1], sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[,
    1], sum))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  musig <- (S * ewu * giepsi/(ssqx2))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig)
  dpmu <- dmusig/pmusig
  sigx1 <- (S * ewu * giepsi/(ssqx2) - dpmu)
  wvsq <- (1 - ewv/sigmasq)
  sigx2 <- (wvsq * ewu/sigmastar)
  sigx3 <- (sigmastar - 0.5 * (ewu * ewv/(ssqx2)))
  wu2sq <- (1 - ewu * gisq/sigmasq)
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (1/(ssqx2) - sigx4 * ewu/ssq)
  pmustar <- pmusig * sigmastar
  sps <- (sigmasq * pmustar)
  sqspmu <- (ssq * pmusig)
  sigx6 <- (ssq * ssqx2)
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = 1/(ssqx2), FUN = "*") -
    sweep(Zigisq, MARGIN = 1, STATS = ewu * sigx3 * giepsi/ssq,
      FUN = "*")
  sigx8 <- (0.5 * sigx2 + sigmastar)
  gradll <- cbind(-(0.5 * (sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv,
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = 2 * (S^2 *
    ewu * giepsi/sigmasq)/ewv, FUN = "*")) + sweep(Xgi, MARGIN = 1,
    STATS = S * dmusig * ewu/sps, FUN = "*")), sweep(uHvar,
    MARGIN = 1, STATS = (0.5 * wu2sq + S * sigx5 * ewu *
      sigx1 * giepsi - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = (0.5 * wvsq + S * sigx8 * dmusig * ewu * ewv *
      giepsi/sqspmu - (0.5 * (2 * (S^2 * sigx8 * ewu^2 *
      ewv * giepsi^2/sigx6) - epsilon_isq/ewv) + 0.5 *
      TT)), FUN = "*"), sweep(sigZ1, MARGIN = 1, STATS = S *
    sigx1 * ewu, FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = 0.5/sigmasq *
    ewu, FUN = "*"))
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phesshalfnormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, ngZGvar, gHvar,
  wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + ngZGvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[,
    1], sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xzigi <- list()
  for (i in 1:ngZGvar) {
    Xzigi[[i]] <- apply(sweep(-Xvar, MARGIN = 1, STATS = gHvar[,
      i] * git, FUN = "*"), 2, function(x) tapply(x, pindex[,
      1], sum))
  }
  Zisqgiepsi <- list()
  for (i in 1:ngZGvar) {
    Zisqgiepsi[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = gHvar[,
      i] * git * epsilon_it, FUN = "*"), 2, function(x) tapply(x,
      pindex[, 1], sum))
  }
  Zisqgisq <- list()
  for (i in 1:ngZGvar) {
    Zisqgisq[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = 4 *
      gHvar[, i] * git^2, FUN = "*"), 2, function(x) tapply(x,
      pindex[, 1], sum))
  }
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) tapply(x, pindex[,
      1], sum))
  }
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  ssqx3 <- sigmasq^2 * sigmastar
  musig <- (S * ewu * giepsi/(ssqx2))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig)
  dpmu <- dmusig/pmusig
  sigx1 <- (S * ewu * giepsi/(ssqx2) - dpmu)
  wvsq <- (1 - ewv/sigmasq)
  sigx2 <- (wvsq * ewu/sigmastar)
  sigx3 <- (sigmastar - 0.5 * (ewu * ewv/(ssqx2)))
  wu2sq <- (1 - ewu * gisq/sigmasq)
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (1/(ssqx2) - sigx4 * ewu/ssq)
  pmustar <- pmusig * sigmastar
  sps <- (sigmasq * pmustar)
  sqspmu <- (ssq * pmusig)
  sigx6 <- (ssq * ssqx2)
  sigx8 <- (0.5 * sigx2 + sigmastar)
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = 1/(ssqx2), FUN = "*") -
    sweep(Zigisq, MARGIN = 1, STATS = ewu * sigx3 * giepsi/ssq,
      FUN = "*")
  sigx9 <- (0.5 + 0.5 * (ewv/sigmasq) - 0.5 * (0.5 * wvsq +
    ewv/sigmasq))
  sigx10 <- (ewu * gisq/sigmasq)
  sigx11 <- ((ewu * gisq/sigmasq - 1) * ewv/sigmasq + 1 - 0.5 *
    (wu2sq * wvsq))
  sigx12 <- (0.5 * wvsq + ewv/sigmasq)
  ewvsqx2 <- (ewv/sigmasq)
  sigx13 <- (0.5 * sigx12 - 0.5 * ewvsqx2)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + ngZGvar,
    ncol = nXvar + nuZUvar + nvZVvar + ngZGvar)
  hessll[1:nXvar, 1:nXvar] <- sapply(1:nXvar, function(x) {
    crossprod(Xsq[[x]], as.matrix(-wHvar * 0.5/ewv))
  }) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * (S^2 *
    ewu/sigmasq)/ewv, FUN = "*"), Xgi) + crossprod(sweep(Xgi,
    MARGIN = 1, STATS = -wHvar * (S^2 * dmusig * (dmusig/sps^2 -
      S * giepsi/(sigmasq^2 * ewv * pmustar)) * ewu^2),
    FUN = "*"), Xgi)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * sigx5 * ewu * (S * (2/sigmastar -
      dmusig * (dmusig/(pmustar) - S * giepsi/ewv)/pmusig) *
      ewu * giepsi/sigmasq - dpmu), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S * ((ewv - S^2 * ewu * giepsi^2/sigmasq)/sqspmu + S *
    ssq * dmusig * ewu * ewv * giepsi/(sqspmu^2 * ssqx2)) *
    sigx8 * dmusig * ewu, FUN = "*"), vHvar) - crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * 2 * (S^2 * sigx8 * ewu^2 *
      ewv * giepsi/sigx6), FUN = "*"), vHvar) + crossprod(sweep(Xepsi_i,
    MARGIN = 1, STATS = wHvar * 0.5/ewv, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- sapply(1:ngZGvar, function(x) {
    crossprod(Xzigi[[x]], as.matrix(wHvar * S * ewu * (sigx1 *
      1/(ssqx2))))
  }) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * S *
    ewu * (sigx1 * (ewu * sigx3 * 1/ssq)), FUN = "*"), Zigisq) +
    crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * S *
      ewu * (S * (1/sigmastar - dmusig * (dmusig/(pmustar) -
      S * giepsi/ewv)/pmusig) * ewu/sigmasq), FUN = "*"),
      sigZ1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ewu * (S * ((1/(ssqx2) - ((0.5 * sigx10 + 1.5 - 0.5 *
    (0.5 * wu2sq + ewu * gisq/sigmasq)) * wu2sq * ewv/sigmastar +
    (3 * gisq - 2 * (sigx4^2 * ewu * sigmasq/ssq)) * sigmastar) *
    ewu/ssq) * sigx1 + S * (1 - dmusig * (dpmu - S * ewu *
    giepsi/(ssqx2))/pmusig) * sigx5^2 * ewu * giepsi) * giepsi -
    0.5 * (wu2sq * gisq/sigmasq)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (0.5 * (gisq/sigmasq^2) +
      S * (((((0.5 * (wu2sq * ewv) - S^2 * sigx8 * sigx5 *
        ewu * giepsi^2)/sigmasq + 0.5 * sigx11 + 0.5 *
        wvsq) * ewu/sigmastar + sigmastar)/sqspmu - sigx8 *
        (2 * (sigx4 * sigmasq * pmustar) - S * ssq *
          sigx5 * dmusig * giepsi) * ewu/sqspmu^2) *
        dmusig - S * (((0.5 * (wu2sq * ewv/sigmasq) +
        0.5 * sigx11) * ewu/sigmastar + 2 * sigx8)/sigx6 -
        ((ssq * gisq + 2 * (sigx4 * ssqx3)) * sigmastar +
          0.5 * (ssq * wu2sq * ewv/sigmastar)) * sigx8 *
          ewu/sigx6^2) * ewu * giepsi) * giepsi) * ewu *
      ewv, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ewu * S * sigx5 * sigx1,
    FUN = "*"), Zigiepsi) + crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ewu * S * ewu * S * (1 - dmusig * (dpmu -
      S * ewu * giepsi/(ssqx2))/pmusig) * sigx5 * giepsi,
    FUN = "*"), sigZ1) - crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ewu * S * ewu * (2 * sigmastar - ((0.25 *
      wu2sq + 0.5 + 0.5 * sigx10) * ewv/(ssqx2) + 2 * (sigx4 *
      ssqx2 * sigx3/ssq)) * ewu) * sigx1/ssq * giepsi,
    FUN = "*"), Zigisq) - crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ewu * 0.5 * (wu2sq/sigmasq), FUN = "*"),
    Zigisq)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (ewv * (S * (((sigx9 * wvsq +
      S^2 * sigx8^2 * ewu * ewv * giepsi^2/(ssq * sigmasq)) *
      ewu/sigmastar + sigmastar)/sqspmu - sigx8^2 * (2 *
      sps + S * dmusig * ewu * giepsi) * ewv/sqspmu^2) *
      dmusig * ewu * giepsi - 0.5 * (wvsq/sigmasq)) - 0.5 *
      (2 * (S^2 * ((sigx9 * wvsq * ewu/sigmastar + sigmastar)/sigx6 -
        ((ssq + 2 * (sigx8 * ssqx3)) * sigmastar + 0.5 *
          (ssq * wvsq * ewu/sigmastar)) * sigx8 * ewv/sigx6^2) *
        ewu^2 * ewv * giepsi^2) + epsilon_isq/ewv)),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      ngZGvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    0.5/sigmasq^2 * ewu * ewv, FUN = "*"), Zigisq) - S *
    (crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((0.5 *
      (wvsq/(ssqx2)) - 0.5 * (1/(ssqx2) - sigx8 * ewv/ssq)) *
      ewu - 2 * (sigx8 * ssqx2 * sigx3/ssq)) * ewu * giepsi *
      sigx1/ssq * ewu * ewv, FUN = "*"), Zigisq) + crossprod(sweep(vHvar,
      MARGIN = 1, STATS = wHvar * sigx8 * sigx1/ssq * ewu *
        ewv, FUN = "*"), Zigiepsi) + crossprod(sweep(vHvar,
      MARGIN = 1, STATS = wHvar * S * sigx8 * (1 + dmusig *
        sigx1/pmusig) * ewu * giepsi/ssq * ewu * ewv,
      FUN = "*"), sigZ1))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + ngZGvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- S * (sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgiepsi[[x]], as.matrix(wHvar * ewu *
        sigx1/(ssqx2)))
    }) - (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
    ewu * sigx1 * (0.5 * (sigx3/ssq) - 0.5/(ssqx3)) * ewu^2 *
    ewv * giepsi * ewu/ssq, FUN = "*"), Zigisq) + sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar * ewu *
        sigx1 * sigx3 * giepsi * ewu/ssq))
    }) + (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
    ewu * sigx1 * sigx3 * ewu/ssq, FUN = "*"), Zigiepsi) +
    crossprod(sweep(Zigiepsi, MARGIN = 1, STATS = wHvar *
      ewu * sigx1 * sigx3 * ewu/ssq, FUN = "*"), Zigisq) -
    crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * ewu *
      sigx1 * 2 * (ewu * ssqx2 * sigx3 * giepsi/ssq) *
      sigx3 * ewu/ssq, FUN = "*"), Zigisq))) + crossprod(sweep(sigZ1,
    MARGIN = 1, STATS = wHvar * ewu * S * (1 - dmusig * (dpmu -
      S * ewu * giepsi/(ssqx2))/pmusig) * ewu, FUN = "*"),
    sigZ1)) - 0.5 * (sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgisq[[x]], as.matrix(wHvar * ewu/sigmasq))
  }) - crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
    ewu * ewu/sigmasq/sigmasq, FUN = "*"), Zigisq))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
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
#' @param modelType specification of inefficiency model G(t)u_i
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
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
halfnormAlgOpt_gzit <- function(start, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar,
  modelType, gHvar, ngZGvar, Yvar, Xvar, pindex, TT, wHvar_c,
  wHvar_p, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, whichStart, initIter, initAlg, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psthalfnorm_gzit(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      gHvar = gHvar, ngZGvar = ngZGvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      wHvar = wHvar_c, modelType = modelType, ngZGvar = ngZGvar,
      gHvar = gHvar, tol = tol, printInfo = printInfo)
    InitHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(phalfnormlike_gzit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    ngZGvar = ngZGvar, gHvar = gHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel BC92-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = phalfnormlike_gzit,
    grad = pgradhalfnormlike_gzit, hess = phesshalfnormlike_gzit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesshalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesshalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(phalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradhalfnormlike_gzit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phesshalfnormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradhalfnormlike_gzit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
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
      mleObj$hessian <- phesshalfnormlike_gzit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesshalfnormlike_gzit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- phalfnormlike_gzit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradhalfnormlike_gzit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, InitHalf = InitHalf))
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
phalfnormeff_gzit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  eta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + object$nvZVvar +
    object$ngZGvar)]
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
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -exp(Wu) * object$S * giepsi/(exp(Wv) + gisq *
    exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  m <- ifelse(mustar > 0, mustar, 0)
  res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB,
    uUB = uUB, m = m, mustar = mustar, sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  res$m <- res$m * git
  res$uLB <- res$uLB * git
  res$uUB <- res$uUB * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teMO <- exp(-res$m)
    res$teBC <- exp(1/2 * res$sigmastar^2 * git^2 - res$mustar *
      git) * pnorm(res$mustar/res$sigmastar - res$sigmastar *
      git)/pnorm(res$mustar/res$sigmastar)
    res$teBCLB <- exp(-res$uUB)
    res$teBCUB <- exp(-res$uLB)
    res$teBC_reciprocal <- exp(1/2 * res$sigmastar^2 * git^2 +
      res$mustar * git) * pnorm(res$mustar/res$sigmastar +
      res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
