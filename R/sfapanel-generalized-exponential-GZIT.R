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
# Convolution: generalized exponential - normal                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for generalized exponential-normal distribution
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
pgenexponormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, ngZGvar,
  gHvar) {
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
  mustar1 <- -(exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  mustar2 <- -(2 * exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  ll <- log(2) - 1/2 * Wu - (TT - 1)/2 * Wv - 1/2 * log(gisq) -
    (TT - 1)/2 * log(2 * pi) + log(exp(-1/2 * (epsilon_isq/exp(Wv) -
    (mustar1/sigmastar)^2)) * pnorm(mustar1/sigmastar) -
    exp(-1/2 * (epsilon_isq/exp(Wv) - (mustar2/sigmastar)^2)) *
      pnorm(mustar2/sigmastar))
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for generalized exponential-normal distribution
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
pstgenexponorm_gzit <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, ngZGvar, gHvar, S, wHvar,
  modelType, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstgenexponorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGenExpo <- NULL
  } else {
    cat("Initialization: SFA + generalized-exponential-normal distribution...\n")
    initGenExpo <- maxLik::maxLik(logLik = cgenexponormlike,
      start = cstgenexponorm(olsObj = olsObj, epsiRes = epsiRes,
        S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
        nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]),
      grad = cgradgenexponormlike, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initGenExpo$estimate
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
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Gradient of the likelihood function ----------
#' gradient for generalized exponential-normal distribution
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
pgradgenexponormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, ngZGvar,
  gHvar) {
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
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmastar <- (sqrt(ewv/gisq) * gisq)
  mustar1 <- (2 * (ewv/ewu_h) + S * giepsi)
  musig1 <- (mustar1/sigmastar)
  pmusig1 <- pnorm(-musig1)
  dmusig1 <- dnorm(-musig1, 0, 1)
  mustar2 <- (ewv/ewu_h + S * giepsi)
  musig2 <- (mustar2/sigmastar)
  pmusig2 <- pnorm(-musig2)
  dmusig2 <- dnorm(-musig2, 0, 1)
  expo1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  expo2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx1 <- (dmusig2 * ewv/sqrt(ewv/gisq))
  sigx2 <- (0.5 * sigx1 - 0.5 * (mustar2 * pmusig2))
  sigx3 <- (dmusig1 * ewv/sqrt(ewv/gisq) - mustar1 * pmusig1)
  sigx4 <- ((expo2 * pmusig2 - expo1 * pmusig1) * ewu_h * gisq)
  sigx5 <- (sigx2 * expo2 - sigx3 * expo1)
  sigx6 <- (1/(ewu_h * gisq) - 0.5 * (mustar2/sigmastar^2))
  sigx7 <- (2 * (sigx6 * mustar2) + epsilon_isq/ewv)
  sigx8 <- (0.5 * (sigx7 * pmusig2) - sigx6 * dmusig2 * ewv/sqrt(ewv/gisq))
  sigx11 <- (2/(ewu_h * gisq) - 0.5 * (mustar1/sigmastar^2))
  sigx9 <- (mustar1 * sigx11)
  sigx10 <- (2 * sigx9 + epsilon_isq/ewv)
  sigx12 <- (0.5 * (sigx10 * pmusig1) - sigx11 * dmusig1 *
    ewv/sqrt(ewv/gisq))
  sigx13 <- (sigx8 * expo2 - sigx12 * expo1)
  sigx14 <- (expo2 * pmusig2 - expo1 * pmusig1)
  sigx15 <- (mustar2 * pmusig2/sigmastar - dmusig2)
  sigx16 <- (sqrt(ewv/gisq) - 0.5 * (ewv/sigmastar))
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar,
    FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = mustar2 *
    sigx16/sigmastar^2, FUN = "*")
  sigx18 <- (mustar1 * pmusig1/sigmastar - dmusig1)
  sigZ2 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar,
    FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = mustar1 *
    sigx16/sigmastar^2, FUN = "*")
  sigZ3 <- sweep(sigZ1, MARGIN = 1, STATS = sigx15 * expo2,
    FUN = "*") - sweep(sigZ2, MARGIN = 1, STATS = sigx18 *
    expo1, FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar1/gisq),
    FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar2/gisq),
    FUN = "*")
  Xsig3 <- sweep(Xepsi_i - 2 * Xsig1, MARGIN = 1, STATS = (pmusig1/ewv),
    FUN = "*")
  Xsig4 <- sweep((Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = (pmusig2/ewv),
    FUN = "*")
  Xsig5 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sigmastar,
    FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sigmastar,
    FUN = "*")
  gradll <- cbind(sweep((0.5 * Xsig3 + Xsig5), MARGIN = 1,
    STATS = expo1/sigx14, FUN = "*") - sweep((0.5 * Xsig4 +
    Xsig6), MARGIN = 1, STATS = expo2/sigx14, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx5/sigx4 - 0.5),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx13/sigx14 -
      0.5 * (TT - 1)), FUN = "*"), sweep(sigZ3, MARGIN = 1,
      STATS = 1/sigx14, FUN = "*") - sweep(Zigisq, MARGIN = 1,
      STATS = 0.5/gisq, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for generalized exponential-normal distribution
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
phessgenexponormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar, ngZGvar,
  gHvar) {
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
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) tapply(x, pindex[,
      1], sum))
  }
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmastar <- (sqrt(ewv/gisq) * gisq)
  mustar1 <- (2 * (ewv/ewu_h) + S * giepsi)
  musig1 <- (mustar1/sigmastar)
  pmusig1 <- pnorm(-musig1)
  dmusig1 <- dnorm(-musig1, 0, 1)
  mustar2 <- (ewv/ewu_h + S * giepsi)
  musig2 <- (mustar2/sigmastar)
  pmusig2 <- pnorm(-musig2)
  dmusig2 <- dnorm(-musig2, 0, 1)
  expo1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  expo2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx1 <- (dmusig2 * ewv/sqrt(ewv/gisq))
  sigx2 <- (0.5 * sigx1 - 0.5 * (mustar2 * pmusig2))
  sigx3 <- (dmusig1 * ewv/sqrt(ewv/gisq) - mustar1 * pmusig1)
  sigx4 <- ((expo2 * pmusig2 - expo1 * pmusig1) * ewu_h * gisq)
  sigx5 <- (sigx2 * expo2 - sigx3 * expo1)
  sigx6 <- (1/(ewu_h * gisq) - 0.5 * (mustar2/sigmastar^2))
  sigx7 <- (2 * (sigx6 * mustar2) + epsilon_isq/ewv)
  sigx8 <- (0.5 * (sigx7 * pmusig2) - sigx6 * dmusig2 * ewv/sqrt(ewv/gisq))
  sigx11 <- (2/(ewu_h * gisq) - 0.5 * (mustar1/sigmastar^2))
  sigx9 <- (mustar1 * sigx11)
  sigx10 <- (2 * sigx9 + epsilon_isq/ewv)
  sigx12 <- (0.5 * (sigx10 * pmusig1) - sigx11 * dmusig1 *
    ewv/sqrt(ewv/gisq))
  sigx13 <- (sigx8 * expo2 - sigx12 * expo1)
  sigx14 <- (expo2 * pmusig2 - expo1 * pmusig1)
  sigx15 <- (mustar2 * pmusig2/sigmastar - dmusig2)
  sigx16 <- (sqrt(ewv/gisq) - 0.5 * (ewv/sigmastar))
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar,
    FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = mustar2 *
    sigx16/sigmastar^2, FUN = "*")
  sigx18 <- (mustar1 * pmusig1/sigmastar - dmusig1)
  sigZ2 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar,
    FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = mustar1 *
    sigx16/sigmastar^2, FUN = "*")
  sigZ3 <- sweep(sigZ1, MARGIN = 1, STATS = sigx15 * expo2,
    FUN = "*") - sweep(sigZ2, MARGIN = 1, STATS = sigx18 *
    expo1, FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar1/gisq),
    FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar2/gisq),
    FUN = "*")
  Xsig3 <- sweep(Xepsi_i - 2 * Xsig1, MARGIN = 1, STATS = (pmusig1/ewv),
    FUN = "*")
  Xsig4 <- sweep((Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = (pmusig2/ewv),
    FUN = "*")
  Xsig5 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sigmastar,
    FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sigmastar,
    FUN = "*")
  dms <- (dmusig2 * mustar2/sigmastar)
  sigx21 <- (pmusig2 - dms)
  sigx22 <- (0.5 * dms - 0.5 * pmusig2)
  sigx23 <- (ewu_h * sqrt(ewv/gisq) * gisq)
  sigx24 <- (ewv/(ewu_h * sigmastar^2))
  sigx25 <- (ewu_h * gisq/(ewu_h * gisq)^2)
  sigx26 <- (ewv/sigmastar^2)
  sigZ4 <- sweep(Zigiepsi, MARGIN = 1, STATS = S * pmusig2,
    FUN = "*") - sweep(sigZ1, MARGIN = 1, STATS = dmusig2 *
    mustar2, FUN = "*")
  sigx29 <- (0.5 * (sigx16/sigmastar^2) - 0.5/(sqrt(ewv/gisq) *
    gisq^2))
  Xsig7 <- sweep((0.5 * Xsig3 + Xsig5), MARGIN = 1, STATS = expo1/sigx14,
    FUN = "*") - sweep((0.5 * Xsig4 + Xsig6), MARGIN = 1,
    STATS = expo2/sigx14, FUN = "*")
  ZZ1 <- list()
  for (i in 1:ngZGvar) {
    ZZ1[[i]] <- sweep(Xzigi[[i]], MARGIN = 1, STATS = 1/sigmastar,
      FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = sigx16 *
      Zigisq[, i]/sigmastar^2, FUN = "*")
  }
  Xsig9 <- sweep((0.5 * Xsig3 + Xsig5), MARGIN = 1, STATS = expo1,
    FUN = "*") - sweep((0.5 * Xsig4 + Xsig6), MARGIN = 1,
    STATS = expo2, FUN = "*")
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + ngZGvar,
    ncol = nXvar + nuZUvar + nvZVvar + ngZGvar)
  hessll[1:nXvar, 1:nXvar] <- 0.5 * (sapply(1:nXvar, function(x) crossprod(Xsq[[x]],
    as.matrix(wHvar * 2 * pmusig1 * expo1/ewv/sigx14))) -
    crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * 2 *
      (S^2/gisq) * pmusig1 * expo1/ewv/sigx14, FUN = "*"),
      Xgi) - crossprod(sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1,
    STATS = wHvar * S * dmusig1/sigmastar * expo1/ewv/sigx14,
    FUN = "*"), Xgi)) - (0.5 * crossprod(sweep((0.5 * Xsig3 +
    Xsig5), MARGIN = 1, STATS = wHvar * expo1/ewv/sigx14,
    FUN = "*"), (Xepsi_i - 2 * Xsig1)) + crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S^2 * mustar1 * dmusig1/(sqrt(ewv/gisq) *
      gisq^2) * expo1/ewv/sigx14, FUN = "*"), Xgi)) - (0.5 *
    (sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar *
      2 * pmusig2 * expo2/ewv/sigx14))) - crossprod(sweep(Xgi,
      MARGIN = 1, STATS = wHvar * 2 * (S^2/gisq) * pmusig2 *
        expo2/ewv/sigx14, FUN = "*"), Xgi) - crossprod(sweep((Xepsi_i -
      2 * Xsig2), MARGIN = 1, STATS = wHvar * S * dmusig2/sigmastar *
      expo2/ewv/sigx14, FUN = "*"), Xgi)) - (0.5 * crossprod(sweep((0.5 *
    Xsig4 + Xsig6), MARGIN = 1, STATS = wHvar * expo2/ewv/sigx14,
    FUN = "*"), (Xepsi_i - 2 * Xsig2)) + crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S^2 * mustar2 * dmusig2/(sqrt(ewv/gisq) *
      gisq^2) * expo2/ewv/sigx14, FUN = "*"), Xgi))) -
    crossprod(sweep(Xsig9, MARGIN = 1, STATS = wHvar/sigx14^2,
      FUN = "*"), Xsig9)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep((Xepsi_i -
    2 * Xsig1), MARGIN = 1, STATS = wHvar * 0.5 * (sigx3/ewv) *
    expo1/sigx4, FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S * pmusig1 * expo1/sigx4, FUN = "*") - (sweep((Xepsi_i -
    2 * Xsig2), MARGIN = 1, STATS = wHvar * 0.5 * (sigx2/ewv) *
    expo2/sigx4, FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S * (0.5 * dms + 0.5 * sigx21) * expo2/sigx4, FUN = "*")) -
    (sweep((0.5 * Xsig3 + Xsig5), MARGIN = 1, STATS = wHvar *
      expo1 * sigx5 * ewu_h * gisq/sigx4^2, FUN = "*") -
      sweep((0.5 * Xsig4 + Xsig6), MARGIN = 1, STATS = wHvar *
        expo2 * sigx5 * ewu_h * gisq/sigx4^2, FUN = "*")),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(0.5 * (sweep(Xgi, MARGIN = 1,
    STATS = wHvar * 2 * (S * (1/(ewu_h * gisq) - mustar2/sigmastar^2)) *
      pmusig2 * expo2, FUN = "*") + sweep(Xepsi_i, MARGIN = 1,
    STATS = wHvar * pmusig2/ewv * expo2, FUN = "*") - sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * sigx7 * dmusig2/sigmastar *
      expo2, FUN = "*")) + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S * (sigx6 * mustar2/(ewv * gisq) + 0.5/sigmastar^2) *
    dmusig2 * ewv/sqrt(ewv/gisq) * expo2, FUN = "*") - sweep((Xepsi_i -
    2 * Xsig2), MARGIN = 1, STATS = wHvar * 0.5 * (sigx8/ewv) *
    expo2, FUN = "*") - (sweep(Xsig7, MARGIN = 1, STATS = wHvar *
    sigx13, FUN = "*") + 0.5 * (sweep(Xgi, MARGIN = 1, STATS = wHvar *
    2 * (S * (2/(ewu_h * gisq) - mustar1/sigmastar^2)) *
    pmusig1 * expo1, FUN = "*") + sweep(Xepsi_i, MARGIN = 1,
    STATS = wHvar * pmusig1/ewv * expo1, FUN = "*") - sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * sigx10 * dmusig1/sigmastar *
      expo1, FUN = "*")) + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S * (mustar1 * sigx11/(ewv * gisq) + 0.5/sigmastar^2) *
    dmusig1 * ewv/sqrt(ewv/gisq) * expo1, FUN = "*") - sweep((Xepsi_i -
    2 * Xsig1), MARGIN = 1, STATS = wHvar * 0.5 * (sigx12/ewv) *
    expo1, FUN = "*")), MARGIN = 1, STATS = 1/sigx14, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * (sigx21/sqrt(ewv/gisq) +
      dmusig2 * mustar2/ewv)/gisq * expo2/sigx14, FUN = "*"),
    sigZ1) - crossprod(sweep((Xepsi_i - 2 * Xsig2), MARGIN = 1,
    STATS = wHvar * 0.5 * (sigx15/ewv) * expo2/sigx14, FUN = "*"),
    sigZ1) + sapply(1:ngZGvar, function(x) crossprod(ZZ1[[x]],
    as.matrix(wHvar * S * sigx15 * expo2/sigx14))) - (crossprod(sweep(Xsig7,
    MARGIN = 1, STATS = wHvar/sigx14, FUN = "*"), sigZ3) +
    crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * S *
      (mustar1 * dmusig1/ewv + (pmusig1 - mustar1 * dmusig1/sigmastar)/sqrt(ewv/gisq))/gisq *
      expo1/sigx14, FUN = "*"), sigZ2) - crossprod(sweep((Xepsi_i -
    2 * Xsig1), MARGIN = 1, STATS = wHvar * 0.5 * (sigx18/ewv) *
    expo1/sigx14, FUN = "*"), sigZ2) + sapply(1:ngZGvar,
    function(x) crossprod(ZZ1[[x]], as.matrix(wHvar * S *
      sigx18 * expo1/sigx14))))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((0.25 * dms - 0.5 * sigx22) * ewv - 0.5 * (sigx2 *
      mustar2/gisq)) * expo2 - expo1 * (ewv * pmusig1 -
      mustar1 * sigx3/gisq))/(sigx14 * ewu_h^2 * gisq) -
      (sigx5/gisq + 0.5 * (sigx14 * ewu_h)) * sigx5 * gisq/sigx4^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (0.5 * (sigx7 * dmusig2 *
      ewv/sigx23) + 2 * (((0.25 * sigx24 - 0.5 * sigx25) *
      mustar2 - 0.5 * (sigx6 * ewv/ewu_h)) * pmusig2)) -
      (((0.25 * sigx26 + 0.5 * (sigx6 * mustar2/gisq))/ewu_h -
        0.5 * sigx25) * dmusig2 * ewv/sqrt(ewv/gisq) +
        0.5 * (sigx8 * mustar2/(ewu_h * gisq)))) * expo2 -
      (sigx13 * sigx5/sigx4 + (0.5 * (sigx10 * dmusig1 *
        ewv/sigx23 + 2 * (((0.5 * sigx24 - sigx25) *
        mustar1 - sigx11 * ewv/ewu_h) * pmusig1)) - (((mustar1 *
        sigx11/gisq + 0.5 * sigx26)/ewu_h - sigx25) *
        dmusig1 * ewv/sqrt(ewv/gisq) + sigx12 * mustar1/(ewu_h *
        gisq))) * expo1))/sigx14, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx22 * ewv/sqrt(ewv/gisq) -
      (0.5 * sigx15 + 0.5 * dmusig2) * mustar2)/gisq *
      expo2/(sigx14 * ewu_h), FUN = "*"), sigZ1) + crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 0.5 * (sigx15 * ewv * sigx16/sigmastar^2) *
      expo2/(sigx14 * ewu_h), FUN = "*"), Zigisq) - (crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((mustar1 * dmusig1/sigmastar -
      pmusig1) * ewv - mustar1^2 * pmusig1/gisq)/sigmastar *
      expo1/(sigx14 * ewu_h), FUN = "*"), sigZ2) + crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx18 * ewv * sigx16/sigmastar^2 *
      expo1/(sigx14 * ewu_h), FUN = "*"), Zigisq) + crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx5/(sigx14 * gisq)/(sigx14 *
      ewu_h), FUN = "*"), sigZ3))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (sigx8 * sigx7) +
      0.5 * ((2 * ((sigx6/ewu_h - 0.5 * ((1/ewu_h - mustar2 *
        gisq/sigmastar^2) * mustar2/sigmastar^2)) * ewv) -
        epsilon_isq/ewv) * pmusig2 - sigx6 * sigx7 *
        dmusig2 * ewv/sqrt(ewv/gisq)) - (1/(ewu_h * gisq) -
      ((sigx6^2 + 0.5/sigmastar^2) * mustar2 + 0.5 * ((1/ewu_h -
        mustar2 * gisq/sigmastar^2) * ewv/sigmastar^2) +
        0.5 * sigx6)) * dmusig2 * ewv/sqrt(ewv/gisq)) *
      expo2 - (sigx13^2/sigx14 + (0.5 * (sigx12 * sigx10) +
      0.5 * ((2 * ((2 * (sigx11/ewu_h) - 0.5 * (mustar1 *
        (2/ewu_h - mustar1 * gisq/sigmastar^2)/sigmastar^2)) *
        ewv) - epsilon_isq/ewv) * pmusig1 - sigx10 *
        sigx11 * dmusig1 * ewv/sqrt(ewv/gisq)) - (2/(ewu_h *
      gisq) - ((sigx11^2 + 0.5/sigmastar^2) * mustar1 +
      0.5 * ((2/ewu_h - mustar1 * gisq/sigmastar^2) * ewv/sigmastar^2) +
      0.5 * sigx11)) * dmusig1 * ewv/sqrt(ewv/gisq)) *
      expo1))/sigx14, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      ngZGvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((pmusig2/ewu_h - sigx6 * dmusig2 * mustar2/sqrt(ewv/gisq))/gisq -
      0.5 * (mustar2 * pmusig2/sigmastar^2)) * ewv/sqrt(ewv/gisq) +
      sigx6 * dmusig2 * mustar2 + 0.5 * (sigx15 * sigx7)) *
    expo2/sigx14, FUN = "*"), sigZ1) - (crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.5/gisq - 0.5 * (1/gisq -
      0.5 * sigx26))/sqrt(ewv/gisq) - sigx16 * gisq/sigmastar^2) *
      mustar2 + sigx16/ewu_h) * sigx15 * ewv/sigmastar^2 *
      expo2/sigx14, FUN = "*"), Zigisq) + crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 0.5 * (S/sqrt(ewv/gisq)) *
      sigx15 * ewv/sigmastar^2 * expo2/sigx14, FUN = "*"),
    Zigiepsi)) - (crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((2 * (pmusig1/ewu_h) - mustar1 * sigx11 * dmusig1/sqrt(ewv/gisq))/gisq -
      0.5 * (mustar1 * pmusig1/sigmastar^2)) * ewv/sqrt(ewv/gisq) +
      mustar1 * sigx11 * dmusig1 + 0.5 * (sigx18 * sigx10)) *
    expo1/sigx14, FUN = "*"), sigZ2) - (crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.5/gisq - 0.5 * (1/gisq -
      0.5 * sigx26))/sqrt(ewv/gisq) - sigx16 * gisq/sigmastar^2) *
      mustar1 + 2 * (sigx16/ewu_h)) * sigx18 * ewv/sigmastar^2 *
      expo1/sigx14, FUN = "*"), Zigisq) + crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 0.5 * (S/sqrt(ewv/gisq)) *
      sigx18 * ewv/sigmastar^2 * expo1/sigx14, FUN = "*"),
    Zigiepsi)) + crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    sigx13/sigx14/sigx14, FUN = "*"), sigZ3))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + ngZGvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- crossprod(sweep(sigZ1,
    MARGIN = 1, STATS = wHvar * mustar2/(ewv * gisq) * mustar2 *
      pmusig2 * expo2/sigx14, FUN = "*") - sweep(Zigisq,
    MARGIN = 1, STATS = wHvar * sigx16/sigmastar^2 * mustar2 *
      pmusig2 * expo2/sigx14, FUN = "*") + sweep(sigZ4,
    MARGIN = 1, STATS = wHvar * expo2/sigmastar/sigx14, FUN = "*"),
    sigZ1) + (S * (sapply(1:ngZGvar, function(x) crossprod(Zisqgiepsi[[x]],
    as.matrix(wHvar * sigx15/sigmastar * expo2/sigx14))) -
    crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx15 *
      sigx16/sigmastar^2 * expo2/sigx14, FUN = "*"), Zigiepsi)) -
    (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
      sigx15 * (sigx29 * ewv * mustar2)/sigmastar^2 * expo2/sigx14,
      FUN = "*") + sweep(Zigiepsi, MARGIN = 1, STATS = wHvar *
      sigx15 * (S * sigx16)/sigmastar^2 * expo2/sigx14,
      FUN = "*"), Zigisq) + (sapply(1:ngZGvar, function(x) crossprod(Zisqgisq[[x]],
      as.matrix(wHvar * sigx15 * mustar2 * sigx16/sigmastar^2 *
        expo2/sigx14))) - crossprod(sweep(Zigisq, MARGIN = 1,
      STATS = wHvar * sigx15 * mustar2 * sigx16 * 2 * (sqrt(ewv/gisq) *
        sigx16 * gisq/sigmastar^2)/sigmastar^2 * expo2/sigx14,
      FUN = "*"), Zigisq)))) - (crossprod(sweep(sigZ2,
    MARGIN = 1, STATS = wHvar * mustar1/(ewv * gisq) * mustar1 *
      pmusig1 * expo1/sigx14, FUN = "*") - sweep(Zigisq,
    MARGIN = 1, STATS = wHvar * sigx16/sigmastar^2 * mustar1 *
      pmusig1 * expo1/sigx14, FUN = "*") + sweep(Zigiepsi,
    MARGIN = 1, STATS = wHvar * S * pmusig1/sigmastar * expo1/sigx14,
    FUN = "*") - sweep(sigZ2, MARGIN = 1, STATS = wHvar *
    mustar1 * dmusig1/sigmastar * expo1/sigx14, FUN = "*"),
    sigZ2) + (S * (sapply(1:ngZGvar, function(x) crossprod(Zisqgiepsi[[x]],
    as.matrix(wHvar * sigx18/sigmastar * expo1/sigx14))) -
    crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx18 *
      sigx16/sigmastar^2 * expo1/sigx14, FUN = "*"), Zigiepsi)) -
    (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
      sigx18 * sigx29 * mustar1 * ewv/sigmastar^2 * expo1/sigx14,
      FUN = "*") + sweep(Zigiepsi, MARGIN = 1, STATS = wHvar *
      sigx18 * S * sigx16/sigmastar^2 * expo1/sigx14, FUN = "*"),
      Zigisq) + (sapply(1:ngZGvar, function(x) crossprod(Zisqgisq[[x]],
      as.matrix(wHvar * sigx18 * mustar1 * sigx16/sigmastar^2 *
        expo1/sigx14))) - crossprod(sweep(Zigisq, MARGIN = 1,
      STATS = wHvar * sigx18 * mustar1 * sigx16 * 2 * (sqrt(ewv/gisq) *
        sigx16 * gisq/sigmastar^2/sigmastar^2) * expo1/sigx14,
      FUN = "*"), Zigisq)))) + crossprod(sweep(sigZ3, MARGIN = 1,
    STATS = wHvar/sigx14^2, FUN = "*"), sigZ3)) - 0.5 * (sapply(1:ngZGvar,
    function(x) crossprod(Zisqgisq[[x]], as.matrix(wHvar/gisq))) -
    crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar/gisq^2,
      FUN = "*"), Zigisq))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for generalized exponential-normal distribution
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
genexponormAlgOpt_gzit <- function(start, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar,
  gHvar, ngZGvar, Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p,
  modelType, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, whichStart, initIter, initAlg, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstgenexponorm_gzit(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      gHvar = gHvar, ngZGvar = ngZGvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar_c, modelType = modelType, ngZGvar = ngZGvar,
      gHvar = gHvar, whichStart = whichStart, initIter = initIter,
      initAlg = initAlg, tol = tol, printInfo = printInfo)
    InitGenExpo <- start_st$initGenExpo
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(pgenexponormlike_gzit(startVal, nXvar = nXvar,
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
      -sum(pgenexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pgenexponormlike_gzit,
    grad = pgradgenexponormlike_gzit, hess = phessgenexponormlike_gzit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pgenexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pgenexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessgenexponormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) {
      -sum(pgenexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessgenexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(pgenexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradgenexponormlike_gzit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phessgenexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradgenexponormlike_gzit(mleObj$par,
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
      mleObj$hessian <- phessgenexponormlike_gzit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessgenexponormlike_gzit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- pgenexponormlike_gzit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradgenexponormlike_gzit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, InitGenExpo = InitGenExpo))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for generalized exponential-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
pgenexponormeff_gzit <- function(object, level) {
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
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar1 <- -(exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  mustar2 <- -(2 * exp(Wv)/(gisq * exp(Wu/2)) + object$S *
    giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  a <- mustar1/sigmastar
  b <- mustar2/sigmastar
  A <- -1/2 * (epsilon_isq/exp(Wv) - a^2)
  B <- -1/2 * (epsilon_isq/exp(Wv) - b^2)
  u <- (exp(A) * (dnorm(a) * sigmastar + mustar1 * pnorm(a)) -
    exp(B) * (dnorm(b) * sigmastar + mustar2 * pnorm(b)))/(exp(A) *
    pnorm(a) - exp(B) * pnorm(b))
  res <- data.frame(levels(pindex[, 1]), u = u, mustar1 = mustar1,
    mustar2 = mustar2, sigmastar = sigmastar, A = A, B = B,
    a = a, b = b)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teBC <- (exp(res$A) * exp(1/2 * res$sigmastar^2 *
      git^2 - res$mustar1 * git) * pnorm(res$a - res$sigmastar *
      git) - exp(res$B) * exp(1/2 * res$sigmastar^2 * git^2 -
      res$mustar2 * git) * pnorm(res$b - res$sigmastar *
      git))/(exp(res$A) * pnorm(res$a) - exp(res$B) * pnorm(res$b))
    res$teBC_reciprocal <- (exp(res$A) * exp(1/2 * res$sigmastar^2 *
      git^2 + res$mustar1 * git) * pnorm(res$a + res$sigmastar *
      git) - exp(res$B) * exp(1/2 * res$sigmastar^2 * git^2 +
      res$mustar2 * git) * pnorm(res$b + res$sigmastar *
      git))/(exp(res$A) * pnorm(res$a) - exp(res$B) * pnorm(res$b))
  }
  res$mustar1 <- NULL
  res$mustar2 <- NULL
  res$A <- NULL
  res$B <- NULL
  res$a <- NULL
  res$b <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
