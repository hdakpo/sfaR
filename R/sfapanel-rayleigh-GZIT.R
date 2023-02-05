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
# Convolution: rayleigh - normal                                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for rayleigh-normal distribution
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
praynormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  mustar <- -exp(Wu) * S * giepsi/(exp(Wv) + gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  ll <- 1/2 * (mustar/sigmastar)^2 - epsilon_isq/(2 * exp(Wv)) -
    (TT - 1)/2 * log(2 * pi) - (TT - 1)/2 * Wv - 1/2 * Wu -
    1/2 * log(exp(Wv) + gisq * exp(Wu)) + log(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for rayleigh-normal distribution
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
pstraynorm_gzit <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, ngZGvar, gHvar, S, wHvar,
  modelType, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstraynorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initRay <- NULL
  } else {
    cat("Initialization: SFA + rayleigh-normal distribution...\n")
    initRay <- maxLik::maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradraynormlike, method = initAlg,
      control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar)
    Esti <- initRay$estimate
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
  return(list(StartVal = StartVal, initRay = initRay))
}

# Gradient of the likelihood function ----------
#' gradient for rayleigh-normal distribution
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
pgradraynormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  Xgi <- apply(Xgit, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  wuvg <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/wuvg)
  siguv <- (wuvg * sigmastar)
  musig <- (S * ewu * giepsi/siguv)
  dmusig <- dnorm(-musig, 0, 1)
  pmusig <- pnorm(-musig)
  dsv <- (dmusig * ewv/sigmastar)
  sigsq <- (sigmastar/ewv - ewu/siguv)
  uuv <- (1 - ewu * gisq/wuvg)
  wvwg <- (1 - ewv/wuvg)
  sigx1 <- (pmusig + S * dmusig * sigsq * giepsi)
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * giepsi/wuvg)
  sigx3 <- (S * giepsi/ewv - sigx1/sigx2)
  sigx4 <- (0.5 * dsv - S * pmusig * giepsi)
  sigx5 <- (0.5 * (uuv * ewv/sigmastar) + sigmastar * gisq)
  sigx6 <- (1/siguv - sigx5 * ewu/siguv^2)
  sigx7 <- (sigx4 * uuv/sigx2 + S^2 * sigx6 * ewu * giepsi^2/sigmastar -
    0.5 * gisq)
  sigx8 <- (wvwg * dmusig/sigmastar)
  sigx9 <- (0.5 * sigx8 + S * pmusig * giepsi/wuvg)
  sigx10 <- (wvwg * ewu/sigmastar)
  sigsig <- (0.5 * sigx10 + sigmastar)
  sigsig2 <- (siguv^2 * sigmastar)
  sigx11 <- (sigx9/sigx2 - S^2 * sigsig * ewu * giepsi^2/sigsig2)
  sigxu <- (sigx11 * ewu - 0.5)
  sigx12 <- (sigmastar - 0.5 * (ewu * ewv/siguv))
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = 1/siguv, FUN = "*") -
    sweep(Zigisq, MARGIN = 1, STATS = ewu * sigx12 * giepsi/siguv^2,
      FUN = "*")
  sigZ2 <- sweep(Zigisq, MARGIN = 1, STATS = dmusig * ewv/siguv,
    FUN = "*")
  sigZ3 <- 0.5 * sigZ2 + sweep(sigZ1, MARGIN = 1, STATS = S^2 *
    dmusig * giepsi, FUN = "*")
  sigZ4 <- sweep(Zigisq, MARGIN = 1, STATS = pmusig/wuvg, FUN = "*") +
    sweep(sigZ1, MARGIN = 1, STATS = S * dmusig, FUN = "*")
  sigZ5 <- sweep(Zigiepsi, MARGIN = 1, STATS = pmusig, FUN = "*") -
    sweep(sigZ4, MARGIN = 1, STATS = ewu * giepsi, FUN = "*")
  sigZ6 <- sweep(sigZ3, MARGIN = 1, STATS = ewu/sigx2, FUN = "*") +
    sweep(sigZ5, MARGIN = 1, STATS = S/sigx2, FUN = "*") +
    0.5 * Zigisq
  sigZ7 <- sweep(sigZ1, MARGIN = 1, STATS = S^2 * ewu * giepsi/sigmastar,
    FUN = "*") - sigZ6
  sigx13 <- epsilon_isq/(2 * ewv)^2
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * ewu *
    sigx3/wuvg, FUN = "*") - sweep(Xepsi_i, MARGIN = 1, STATS = 1/(2 *
    ewv), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx7 *
    ewu/wuvg - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = ((sigxu/wuvg + 2 * (sigx13)) * ewv - 0.5 * (TT -
      1)), FUN = "*"), sweep(sigZ7, MARGIN = 1, STATS = ewu/wuvg,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for rayleigh-normal distribution
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
phessraynormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  Xgi <- apply(Xgit, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, index(data)[[1]],
    sum))
  Xzigi <- list()
  for (i in 1:ngZGvar) {
    Xzigi[[i]] <- apply(sweep(-Xvar, MARGIN = 1, STATS = gHvar[,
      i] * git, FUN = "*"), 2, function(x) tapply(x, index(data)[[1]],
      sum))
  }
  Zisqgiepsi <- list()
  for (i in 1:ngZGvar) {
    Zisqgiepsi[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = gHvar[,
      i] * git * epsilon_it, FUN = "*"), 2, function(x) tapply(x,
      index(data)[[1]], sum))
  }
  Zisqgisq <- list()
  for (i in 1:ngZGvar) {
    Zisqgisq[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = 4 *
      gHvar[, i] * git^2, FUN = "*"), 2, function(x) tapply(x,
      index(data)[[1]], sum))
  }
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) tapply(x, index(data)[[1]],
      sum))
  }
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  wuvg <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/wuvg)
  siguv <- (wuvg * sigmastar)
  musig <- (S * ewu * giepsi/siguv)
  dmusig <- dnorm(-musig, 0, 1)
  pmusig <- pnorm(-musig)
  dsv <- (dmusig * ewv/sigmastar)
  sigsq <- (sigmastar/ewv - ewu/siguv)
  uuv <- (1 - ewu * gisq/wuvg)
  wvwg <- (1 - ewv/wuvg)
  sigx1 <- (pmusig + S * dmusig * sigsq * giepsi)
  sigx2 <- (dmusig * sigmastar - S * ewu * pmusig * giepsi/wuvg)
  sigx3 <- (S * giepsi/ewv - sigx1/sigx2)
  sigx4 <- (0.5 * dsv - S * pmusig * giepsi)
  sigx5 <- (0.5 * (uuv * ewv/sigmastar) + sigmastar * gisq)
  sigx6 <- (1/siguv - sigx5 * ewu/siguv^2)
  sigx7 <- (sigx4 * uuv/sigx2 + S^2 * sigx6 * ewu * giepsi^2/sigmastar -
    0.5 * gisq)
  sigx8 <- (wvwg * dmusig/sigmastar)
  sigx9 <- (0.5 * sigx8 + S * pmusig * giepsi/wuvg)
  sigx10 <- (wvwg * ewu/sigmastar)
  sigsig <- (0.5 * sigx10 + sigmastar)
  sigsig2 <- (siguv^2 * sigmastar)
  sigx11 <- (sigx9/sigx2 - S^2 * sigsig * ewu * giepsi^2/sigsig2)
  sigxu <- (sigx11 * ewu - 0.5)
  sigx12 <- (sigmastar - 0.5 * (ewu * ewv/siguv))
  sigx13 <- epsilon_isq/(2 * ewv)^2
  sigx14 <- (0.5 * (sigx12/siguv^2) - 0.5/(wuvg^2 * sigmastar))
  sigx15 <- (S * sigsig * dmusig * ewu * giepsi/siguv^2 - pmusig/wuvg)
  sigx16 <- (pmusig * gisq/wuvg + S * sigx6 * dmusig * giepsi)
  sigx17 <- (sigx5 * wuvg * sigmastar * sigx12/siguv^2)
  sigx18 <- ((0.5 * uuv - 0.5)/siguv - 0.5 * sigx6)
  sigx19 <- ((sigx18 * ewv - 2 * sigx17) * ewu + sigmastar)
  sigx20 <- (ewu * gisq/wuvg)
  sigx21 <- (siguv^2 * wuvg * sigmastar)
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = 1/siguv, FUN = "*") -
    sweep(Zigisq, MARGIN = 1, STATS = ewu * sigx12 * giepsi/siguv^2,
      FUN = "*")
  sigZ2 <- sweep(Zigisq, MARGIN = 1, STATS = dmusig * ewv/siguv,
    FUN = "*")
  sigZ3 <- 0.5 * sigZ2 + sweep(sigZ1, MARGIN = 1, STATS = S^2 *
    dmusig * giepsi, FUN = "*")
  sigZ4 <- sweep(Zigisq, MARGIN = 1, STATS = pmusig/wuvg, FUN = "*") +
    sweep(sigZ1, MARGIN = 1, STATS = S * dmusig, FUN = "*")
  sigZ5 <- sweep(Zigiepsi, MARGIN = 1, STATS = pmusig, FUN = "*") -
    sweep(sigZ4, MARGIN = 1, STATS = ewu * giepsi, FUN = "*")
  sigZ6 <- sweep(sigZ3, MARGIN = 1, STATS = ewu/sigx2, FUN = "*") +
    sweep(sigZ5, MARGIN = 1, STATS = S/sigx2, FUN = "*") +
    0.5 * Zigisq
  sigZ7 <- sweep(sigZ1, MARGIN = 1, STATS = S^2 * ewu * giepsi/sigmastar,
    FUN = "*") - sigZ6
  sigZ8 <- sweep(Zigisq, MARGIN = 1, STATS = (ewu * giepsi/wuvg),
    FUN = "*")
  sigZ9 <- sweep(Zigisq, MARGIN = 1, STATS = (ewu * wuvg *
    sigmastar * sigx12 * giepsi/siguv^2), FUN = "*")
  sigZ10 <- sweep(sigZ3, MARGIN = 1, STATS = ewu, FUN = "*") +
    S * sigZ5
  sigZ11 <- sweep(Zigisq, MARGIN = 1, STATS = sigx19 * giepsi,
    FUN = "*") + sweep(Zigiepsi, MARGIN = 1, STATS = sigx5,
    FUN = "*")
  sigZ12 <- sweep(sigZ11, MARGIN = 1, STATS = 1/siguv^2, FUN = "*") +
    sweep(sigZ1, MARGIN = 1, STATS = S^2 * sigx6 * ewu *
      giepsi^2/siguv, FUN = "*")
  sigZ13 <- sweep(sigZ1, MARGIN = 1, STATS = S^2 * ewu^2 *
    giepsi^2/siguv, FUN = "*")
  sigZ14 <- sweep(sigZ10, MARGIN = 1, STATS = sigx1/sigx2,
    FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = pmusig,
    FUN = "*")
  sigZ15 <- sweep(sigZ1, MARGIN = 1, STATS = (1 - S^2 * ewu *
    giepsi^2/(wuvg * ewv)), FUN = "*")
  sigZ16 <- sweep(sigZ1, MARGIN = 1, STATS = ewu * S^2 * giepsi^2/ewv,
    FUN = "*") + sweep(Zigisq, MARGIN = 1, STATS = ewu *
    giepsi/siguv, FUN = "*")
  sigZ17 <- sweep(sigZ16, MARGIN = 1, STATS = 1/wuvg, FUN = "*") -
    sweep((0.5 * sigZ8 + 2 * Zigiepsi), MARGIN = 1, STATS = 1/(sigmastar *
      wuvg), FUN = "*")
  sigZ18 <- sweep(Zigisq, MARGIN = 1, STATS = ewu * sigx12 *
    giepsi/siguv^2, FUN = "*")
  sigZ19 <- sweep(sigZ14, MARGIN = 1, STATS = ewu/wuvg, FUN = "*") +
    sweep((sigZ15 + sigZ17 + sigZ18), MARGIN = 1, STATS = S *
      dmusig * ewu, FUN = "*")
  sigZ20 <- sweep(Zigisq, MARGIN = 1, STATS = ((((0.5 * uuv -
    1.5)/siguv - 0.5 * sigx6) * ewv - 2 * sigx17) * ewu +
    3 * sigmastar) * giepsi, FUN = "*")
  sigZ21 <- sweep(Zigisq, MARGIN = 1, STATS = (dmusig * ewv/sigmastar),
    FUN = "*")
  sigZ22 <- sweep(Zigisq, MARGIN = 1, STATS = (sigx16/wuvg) *
    ewu, FUN = "*") + sweep(sigZ12, MARGIN = 1, STATS = (S *
    dmusig) * ewu, FUN = "*")
  sigZ23 <- sweep(Zigiepsi, MARGIN = 1, STATS = S * (sigx5 *
    ewu/siguv^2 - 1/siguv) * dmusig, FUN = "*")
  sigZ24 <- sweep(Zigisq, MARGIN = 1, STATS = ((sigx5 * dmusig *
    ewv/siguv^2 + S^2 * sigx6 * dmusig * giepsi^2/wuvg)),
    FUN = "*")
  sigZ25 <- sweep(Zigiepsi, MARGIN = 1, STATS = (2/siguv) *
    giepsi/sigmastar, FUN = "*") - (sweep(sigZ20, MARGIN = 1,
    STATS = ewu/siguv^2 * giepsi/sigmastar, FUN = "*") +
    sweep(Zigiepsi, MARGIN = 1, STATS = sigx5 * ewu/siguv^2 *
      giepsi/sigmastar, FUN = "*") + sweep(sigZ1, MARGIN = 1,
    STATS = 0.5 * (uuv) * giepsi/sigmastar, FUN = "*"))
  sigZ26 <- sweep(0.5 * sigZ21, MARGIN = 1, STATS = 1/wuvg,
    FUN = "*") - sweep(sigZ10, MARGIN = 1, STATS = sigx4 *
    uuv/sigx2/wuvg, FUN = "*") + sweep(sigZ22 + sigZ23, MARGIN = 1,
    STATS = S * giepsi, FUN = "*") - sweep(Zigisq, MARGIN = 1,
    STATS = pmusig/wuvg * S * giepsi, FUN = "*")
  sigZ27 <- sweep(sigZ26, MARGIN = 1, STATS = 1/sigx2, FUN = "*") -
    (sweep(0.5 * sigZ24, MARGIN = 1, STATS = ewu/sigx2, FUN = "*") +
      sweep(sigZ12, MARGIN = 1, STATS = S^2 * dmusig *
        giepsi * ewu/sigx2, FUN = "*")) + sweep(sigZ7,
    MARGIN = 1, STATS = gisq/wuvg, FUN = "*")
  sigZ28 <- sweep((S^2 * sigZ25 - sigZ27), MARGIN = 1, STATS = ewu/wuvg,
    FUN = "*") - sweep(sigZ6, MARGIN = 1, STATS = 1/wuvg,
    FUN = "*")
  sigZ29 <- sweep(Zigisq, MARGIN = 1, STATS = ((dmusig + S^2 *
    sigsig * dmusig * ewu^2 * ewv * giepsi^2/sigx21)/siguv -
    sigsig * dmusig * ewv/siguv^2), FUN = "*")
  sigZ30 <- sweep(Zigiepsi, MARGIN = 1, STATS = S * sigsig *
    dmusig/siguv^2 * S * giepsi, FUN = "*") - sweep(Zigisq,
    MARGIN = 1, STATS = sigx15/wuvg * S * giepsi, FUN = "*") -
    sweep(sigZ10, MARGIN = 1, STATS = sigx9/(sigx2 * wuvg),
      FUN = "*")
  sigZ31 <- sweep(Zigisq, MARGIN = 1, STATS = ((0.5 * (wvwg/siguv) -
    0.5 * (1/siguv - sigsig * ewv/siguv^2)) * ewu - 2 * (sigsig *
    wuvg * sigmastar * sigx12/siguv^2)) * ewu * giepsi, FUN = "*")
  sigZ32 <- sweep(sigZ31, MARGIN = 1, STATS = ewv/siguv^2 *
    giepsi/sigmastar, FUN = "*") + sweep(Zigiepsi, MARGIN = 1,
    STATS = sigsig * ewv/siguv^2 * giepsi/sigmastar, FUN = "*") +
    sweep(sigZ1, MARGIN = 1, STATS = 0.5 * (wvwg) * giepsi/sigmastar,
      FUN = "*")
  sigZ33 <- sweep((0.5 * (sigZ29) + sigZ30), MARGIN = 1, STATS = ewv/sigx2 *
    ewu * ewu/wuvg, FUN = "*") + sweep(S^2 * sigZ32, MARGIN = 1,
    STATS = ewu * ewu/wuvg, FUN = "*") + sweep(sigZ7, MARGIN = 1,
    STATS = ewv/wuvg * ewu/wuvg, FUN = "*")
  sigZ34 <- sweep(sigZ1, MARGIN = 1, STATS = 2 * (S * dmusig),
    FUN = "*") + sweep(Zigisq, MARGIN = 1, STATS = pmusig/wuvg,
    FUN = "*")
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + ngZGvar,
    ncol = nXvar + nuZUvar + nvZVvar + ngZGvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S^2 * (1/ewv - (((1 - S^2 * ewu * giepsi^2/(wuvg *
      ewv)) * sigsq - ewu/siguv) * dmusig + ewu * sigx1^2/(sigx2 *
      wuvg))/sigx2) * ewu/wuvg, FUN = "*"), Xgi) - sapply(1:nXvar,
    function(x) crossprod(Xsq[[x]], as.matrix(wHvar/ewv)))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * (((sigx4 * sigx1/sigx2 +
      0.5 * (S * dmusig * giepsi/sigmastar)) * ewu/wuvg -
      pmusig) * uuv/sigx2 + 2 * (S * sigx6 * ewu * giepsi/sigmastar)) *
      ewu/wuvg, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xepsi_i, MARGIN = 1, STATS = wHvar *
    2 * (1/(2 * ewv)^2) * ewv, FUN = "*") + sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * (((sigx9 * sigx1/sigx2 - S * (0.5 *
      (wvwg/ewv) + 1/wuvg) * dmusig * giepsi/sigmastar) *
      ewu + pmusig)/(sigx2 * wuvg) - 2 * (S * sigsig *
      ewu * giepsi/sigsig2)) * ewu/wuvg * ewv, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- sapply(1:ngZGvar, function(x) crossprod(Xzigi[[x]],
    as.matrix(wHvar * S * ewu/wuvg * giepsi * S * ewu/(siguv *
      sigmastar)))) - crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * ewu/wuvg * ewu^2 * sigx12/siguv^2 *
      S * giepsi/sigmastar, FUN = "*"), Zigisq) + crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * ewu/wuvg * S * ewu/sigmastar,
    FUN = "*"), sigZ1) - (crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * ewu/(sigx2 * wuvg), FUN = "*"), sigZ19) +
    sapply(1:ngZGvar, function(x) crossprod(Xzigi[[x]], as.matrix(wHvar *
      S * ewu/wuvg * pmusig/sigx2))))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * dsv + ewu * (S^2 * sigx6 * dmusig * giepsi^2 -
      (sigx4 * uuv/sigx2 + gisq) * sigx4/wuvg) - (0.5 *
      (0.5 * (uuv * dmusig * ewv/sigmastar) + S^2 * sigx6 *
        dmusig * ewu * giepsi^2) + S * pmusig * giepsi)) *
      uuv/sigx2 + ewu * (S^2 * (2/siguv - (((0.5 * sigx20 +
      2 - 0.5 * (0.5 * uuv + ewu * gisq/wuvg)) * uuv *
      ewv/sigmastar + (4 * gisq - 2 * (sigx5^2 * ewu *
      wuvg/siguv^2)) * sigmastar) * ewu/siguv^2 + 0.5 *
      (uuv * sigx6))) * giepsi^2/sigmastar - sigx7 * gisq/wuvg) -
      0.5 * gisq) * ewu/wuvg, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (wvwg * dmusig) +
      0.5 * ((dmusig * ewv * gisq/wuvg - S^2 * wvwg * sigx6 *
        dmusig * ewu * giepsi^2/sigmastar) * ewu/wuvg -
        0.5 * (uuv * wvwg * dmusig)))/sigmastar + (S *
      pmusig * giepsi - (sigx9 * sigx4 * uuv/sigx2 + S *
      sigx16 * giepsi) * ewu)/wuvg)/sigx2 - (sigxu * gisq/wuvg +
      S^2 * (((0.5 * (uuv * ewv/wuvg) + 0.5 * ((ewu * gisq/wuvg -
        1) * ewv/wuvg + 1 - 0.5 * (uuv * wvwg)) + 0.5 *
        wvwg) * ewu/sigmastar + sigmastar)/sigsig2 +
        sigsig * (1/sigsig2 - (0.5 * (siguv^2 * uuv/siguv) +
          2 * (sigx5 * ewu)) * ewu * ewv/sigsig2^2)) *
        ewu * giepsi^2)) * ewu * ewv/wuvg, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ewu, FUN = "*"), sigZ28)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (ewv * (S^2 *
      sigsig * dmusig * ewu^2 * giepsi^2/sigsig2 - dmusig)/wuvg -
      0.5 * (wvwg * dmusig)) + 0.5 * dmusig) * wvwg/sigmastar +
      (ewv * (S * sigx15 * giepsi - sigx9^2 * ewu/sigx2) +
        S * pmusig * giepsi)/wuvg)/sigx2 - S^2 * (sigsig *
      (1/sigsig2 - (0.5 * (siguv^2 * wvwg/siguv) + 2 *
        (sigsig * ewv)) * ewu * ewv/sigsig2^2) + (0.5 *
      (ewv/wuvg) - 0.5 * (0.5 * wvwg + ewv/wuvg)) * wvwg *
      wuvg/(siguv^2 * ewv)) * ewu * giepsi^2) * ewu - (sigxu *
      ewv/wuvg + 0.5))/wuvg + (2 - 16 * (ewv^2/(2 * ewv)^2)) *
      sigx13) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      ngZGvar)] <- -(crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar, FUN = "*"), sigZ33))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + ngZGvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- ((S^2 * (((0.5 * crossprod(sweep(sigZ1,
    MARGIN = 1, STATS = wHvar/wuvg * ewu * giepsi/sigmastar *
      ewu * ewu/wuvg, FUN = "*"), Zigisq) - ((crossprod(sweep(Zigisq,
    MARGIN = 1, STATS = wHvar * sigx14 * ewu^2 * ewv * giepsi/siguv^2 *
      ewu * giepsi/sigmastar * ewu * ewu/wuvg, FUN = "*"),
    Zigisq) + sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgisq[[x]], as.matrix(wHvar * sigx12 * giepsi/siguv^2 *
      ewu * giepsi/sigmastar * ewu * ewu/wuvg))
  })) + crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
    sigx12/siguv^2 * ewu * giepsi/sigmastar * ewu * ewu/wuvg,
    FUN = "*"), (2 * Zigiepsi - 2 * sigZ9)))) + sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgiepsi[[x]], as.matrix(wHvar/siguv *
        giepsi/sigmastar * ewu * ewu/wuvg))
    })) + crossprod(sweep(Zigiepsi, MARGIN = 1, STATS = wHvar/sigmastar *
    ewu * ewu/wuvg, FUN = "*"), sigZ1)) - crossprod(sweep(sigZ7,
    MARGIN = 1, STATS = wHvar/wuvg * ewu * ewu/wuvg, FUN = "*"),
    Zigisq)) - (((crossprod(sweep(sigZ10, MARGIN = 1, STATS = wHvar/(sigx2 *
    wuvg) * ewu/sigx2 * ewu/wuvg, FUN = "*"), sigZ10) + 0.5 *
    (((sapply(1:ngZGvar, function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar * dmusig/siguv *
        ewv * ewu/sigx2 * ewu/wuvg))
    }) - crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
      S^2 * dmusig * ewu^2 * giepsi/siguv/siguv * ewv *
      ewu/sigx2 * ewu/wuvg, FUN = "*"), sigZ1)) - crossprod(sweep(Zigisq,
      MARGIN = 1, STATS = wHvar * dmusig * ewu * sigx12/siguv^2 *
        ewv * ewu/sigx2 * ewu/wuvg, FUN = "*"), Zigisq))) +
    S^2 * (crossprod(sweep((Zigiepsi - sigZ13), MARGIN = 1,
      STATS = wHvar * dmusig * ewu/sigx2 * ewu/wuvg, FUN = "*"),
      sigZ1) + (sapply(1:ngZGvar, function(x) {
      crossprod(Zisqgiepsi[[x]], as.matrix(wHvar/siguv *
        dmusig * giepsi * ewu/sigx2 * ewu/wuvg))
    }) - ((crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
      sigx14 * ewu^2 * ewv * giepsi * ewu/siguv^2 * dmusig *
      giepsi * ewu/sigx2 * ewu/wuvg, FUN = "*"), Zigisq) +
      sapply(1:ngZGvar, function(x) {
        crossprod(Zisqgisq[[x]], as.matrix(wHvar * sigx12 *
          giepsi * ewu/siguv^2 * dmusig * giepsi * ewu/sigx2 *
          ewu/wuvg))
      })) + crossprod(sweep((2 * Zigiepsi - 2 * sigZ9),
      MARGIN = 1, STATS = wHvar * sigx12 * ewu/siguv^2 *
        dmusig * giepsi * ewu/sigx2 * ewu/wuvg, FUN = "*"),
      Zigisq))))) + S * (sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgiepsi[[x]], as.matrix(wHvar * pmusig/sigx2 *
      ewu/wuvg))
  }) - (((sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgisq[[x]], as.matrix(wHvar * pmusig/wuvg *
      giepsi * ewu/sigx2 * ewu/wuvg))
  }) - crossprod(sweep(sigZ4, MARGIN = 1, STATS = wHvar * ewu/wuvg *
    giepsi * ewu/sigx2 * ewu/wuvg, FUN = "*"), Zigisq)) +
    (sapply(1:ngZGvar, function(x) {
      crossprod(Zisqgiepsi[[x]], as.matrix(wHvar/siguv *
        S * dmusig * giepsi * ewu/sigx2 * ewu/wuvg))
    }) - (((crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
      sigx14 * ewu^2 * ewv * giepsi/siguv^2 * ewu * S *
      dmusig * giepsi * ewu/sigx2 * ewu/wuvg, FUN = "*"),
      Zigisq) + sapply(1:ngZGvar, function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar * sigx12 *
        giepsi/siguv^2 * ewu * S * dmusig * giepsi *
        ewu/sigx2 * ewu/wuvg))
    })) + crossprod(sweep((2 * Zigiepsi - 2 * sigZ9), MARGIN = 1,
      STATS = wHvar * sigx12/siguv^2 * ewu * S * dmusig *
        giepsi * ewu/sigx2 * ewu/wuvg, FUN = "*"), Zigisq)) +
      crossprod(sweep(sigZ1, MARGIN = 1, STATS = wHvar *
        S^2 * ewu * giepsi/siguv * ewu * S * dmusig *
        giepsi * ewu/sigx2 * ewu/wuvg, FUN = "*"), sigZ1)))) +
    crossprod(sweep(sigZ34, MARGIN = 1, STATS = wHvar * ewu/sigx2 *
      ewu/wuvg, FUN = "*"), Zigiepsi)))) + 0.5 * sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar * ewu/wuvg))
    })))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for rayleigh-normal distribution
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
raynormAlgOpt_gzit <- function(start, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar,
  gHvar, ngZGvar, modelType, Yvar, Xvar, pindex, TT, wHvar_c,
  wHvar_p, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstraynorm_gzit(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      gHvar = gHvar, ngZGvar = ngZGvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar_c, ngZGvar = ngZGvar, gHvar = gHvar,
      initIter = initIter, initAlg = initAlg, modelType = modelType,
      whichStart = whichStart, tol = tol, printInfo = printInfo)
    InitRay <- start_st$initRay
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(praynormlike_gzit(startVal, nXvar = nXvar,
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
      -sum(praynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = praynormlike_gzit,
    grad = pgradraynormlike_gzit, hess = phessraynormlike_gzit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(praynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(praynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessraynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) {
      -sum(praynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_gzit(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessraynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(praynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradraynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phessraynormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradraynormlike_gzit(mleObj$par,
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
      mleObj$hessian <- phessraynormlike_gzit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessraynormlike_gzit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- praynormlike_gzit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradraynormlike_gzit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, InitRay = InitRay))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for rayleigh-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
praynormeff_gzit <- function(object, level) {
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
  u <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 +
    sigmastar^2) * pnorm(mustar/sigmastar))/(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  m <- ifelse(mustar/2 + sqrt(sigmastar^2 + mustar^2/4) > 0,
    mustar/2 + sqrt(sigmastar^2 + mustar^2/4), 0)
  res <- data.frame(levels(pindex[, 1]), u = u, m = m, mustar = mustar,
    sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  res$m <- res$m * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teMO <- exp(-res$m)
    res$teBC <- exp(-res$mustar * git + res$sigmastar^2 *
      git^2/2) * (res$sigmastar * git * dnorm(res$mustar/res$sigmastar -
      res$sigmastar * git) + (res$mustar * git - res$sigmastar^2 *
      git^2) * pnorm(res$mustar/res$sigmastar - res$sigmastar *
      git))/(res$sigmastar * git * dnorm(res$mustar/res$sigmastar) +
      res$mustar * git * pnorm(res$mustar/res$sigmastar))
    res$teBC_reciprocal <- exp(res$mustar * git + res$sigmastar^2 *
      git^2/2) * (res$sigmastar * git * dnorm(res$mustar/res$sigmastar +
      res$sigmastar * git) + (res$mustar * git + res$sigmastar^2 *
      git^2) * pnorm(res$mustar/res$sigmastar + res$sigmastar *
      git))/(res$sigmastar * git * dnorm(res$mustar/res$sigmastar) +
      res$mustar * git * pnorm(res$mustar/res$sigmastar))
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
