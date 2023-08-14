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
# Convolution: generalized exponential - normal                                #
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
pgenexponormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  mustar1 <- -(exp(Wv)/(TT * exp(Wu/2)) + S * epsilon_i/TT)
  mustar2 <- -(2 * exp(Wv)/(TT * exp(Wu/2)) + S * epsilon_i/TT)
  sigmastar <- sqrt(exp(Wv)/TT)
  ll <- log(2) - 1/2 * Wu - (TT - 1)/2 * Wv - 1/2 * log(TT) - (TT - 1)/2 * log(2 *
    pi) + log(exp(-1/2 * (epsilon_isq/exp(Wv) - (mustar1/sigmastar)^2)) * pnorm(mustar1/sigmastar) -
    exp(-1/2 * (epsilon_isq/exp(Wv) - (mustar2/sigmastar)^2)) * pnorm(mustar2/sigmastar))
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
pstgenexponorm_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- cstgenexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGenExpo <- NULL
  } else {
    cat("Initialization: SFA + generalized-exponential-normal distribution...\n")
    initGenExpo <- maxLik::maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradgenexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initGenExpo$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  })
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
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
pgradgenexponormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  X_iM <- apply(-Xvar, 2, function(x) tapply(x, pindex[, 1], sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  epsi_vu <- (ewv/ewu_h + S * epsilon_i)
  wvtt <- ewv/TT
  ssq_v <- (TT * sqrt(wvtt))
  wutt <- TT * ewu_h
  musig <- (epsi_vu/ssq_v)
  dmusig <- dnorm(-musig)
  pmusig <- pnorm(-musig)
  epsi_vu2 <- (2 * (ewv/ewu_h) + S * epsilon_i)
  musig2 <- (epsi_vu2/ssq_v)
  pmusig2 <- pnorm(-musig2)
  dmusig2 <- dnorm(-musig2, 0, 1)
  sigx1_1 <- (epsi_vu/ssq_v^2)
  sigx1_2 <- (epsi_vu2/ssq_v^2)
  sigx2_1 <- (1/(wutt) - 0.5 * sigx1_1)
  sigx2_2 <- (epsi_vu2 * (2/(wutt) - 0.5 * sigx1_2))
  sigx3_1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig)^2)))
  sigx3_2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx4 <- (sigx3_1 * pmusig - sigx3_2 * pmusig2)
  ttsq <- (TT * sigx4 * ewu_h)
  epsivusq1 <- (2 * (sigx2_1 * epsi_vu) + epsilon_isq/ewv)
  epsivusq2 <- (2 * sigx2_2 + epsilon_isq/ewv)
  wusig <- (2/(wutt) - 0.5 * sigx1_2)
  dwv1 <- dmusig * ewv/sqrt(wvtt)
  dwv2 <- dmusig2 * ewv/sqrt(wvtt)
  sigx5_1 <- (0.5 * (epsivusq1 * pmusig) - sigx2_1 * dwv1)
  sigx5_2 <- (0.5 * (epsivusq2 * pmusig2) - wusig * dwv2)
  sigx6 <- (sigx5_1 * sigx3_1 - sigx5_2 * sigx3_2)
  dw1pmu <- (0.5 * (dwv1) - 0.5 * (epsi_vu * pmusig))
  dw2pmu <- (dwv2 - epsi_vu2 * pmusig2)
  sigx7 <- (dw1pmu * sigx3_1 - dw2pmu * sigx3_2)
  Xsig1 <- sweep(2 * Xepsi_i, MARGIN = 1, STATS = pmusig2/ewv, FUN = "*") - sweep(X_iM,
    MARGIN = 1, STATS = 2 * (S * epsi_vu2/TT) * pmusig2/ewv, FUN = "*")
  Xsig2 <- sweep(X_iM, MARGIN = 1, STATS = S * dmusig2/ssq_v, FUN = "*")
  Xsig3 <- sweep(2 * Xepsi_i, MARGIN = 1, STATS = pmusig/ewv, FUN = "*") - sweep(X_iM,
    MARGIN = 1, STATS = 2 * (S * epsi_vu/TT) * pmusig/ewv, FUN = "*")
  Xsig4 <- sweep(X_iM, MARGIN = 1, STATS = S * dmusig/ssq_v, FUN = "*")
  Xsig5 <- sweep((0.5 * Xsig1 + Xsig2), MARGIN = 1, STATS = sigx3_2, FUN = "*") -
    sweep((0.5 * Xsig3 + Xsig4), MARGIN = 1, STATS = sigx3_1, FUN = "*")
  gradll <- cbind(sweep(Xsig5, MARGIN = 1, STATS = 1/sigx4, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = (sigx7/ttsq - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = (sigx6/sigx4 - 0.5 * (TT - 1)), FUN = "*"))
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
phessgenexponormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  X_iM <- apply(-Xvar, 2, function(x) tapply(x, pindex[, 1], sum))
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  epsi_vu <- (ewv/ewu_h + S * epsilon_i)
  wvtt <- ewv/TT
  ssq_v <- (TT * sqrt(wvtt))
  wutt <- TT * ewu_h
  musig <- (epsi_vu/ssq_v)
  dmusig <- dnorm(-musig)
  pmusig <- pnorm(-musig)
  epsi_vu2 <- (2 * (ewv/ewu_h) + S * epsilon_i)
  musig2 <- (epsi_vu2/ssq_v)
  pmusig2 <- pnorm(-musig2)
  dmusig2 <- dnorm(-musig2, 0, 1)
  sigx1_1 <- (epsi_vu/ssq_v^2)
  sigx1_2 <- (epsi_vu2/ssq_v^2)
  sigx2_1 <- (1/(wutt) - 0.5 * sigx1_1)
  sigx2_2 <- (epsi_vu2 * (2/(wutt) - 0.5 * sigx1_2))
  sigx3_1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig)^2)))
  sigx3_2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx4 <- (sigx3_1 * pmusig - sigx3_2 * pmusig2)
  ttsq <- (TT * sigx4 * ewu_h)
  epsivusq1 <- (2 * (sigx2_1 * epsi_vu) + epsilon_isq/ewv)
  epsivusq2 <- (2 * sigx2_2 + epsilon_isq/ewv)
  wusig <- (2/(wutt) - 0.5 * sigx1_2)
  dwv1 <- dmusig * ewv/sqrt(wvtt)
  dwv2 <- dmusig2 * ewv/sqrt(wvtt)
  sigx5_1 <- (0.5 * (epsivusq1 * pmusig) - sigx2_1 * dwv1)
  sigx5_2 <- (0.5 * (epsivusq2 * pmusig2) - wusig * dwv2)
  sigx6 <- (sigx5_1 * sigx3_1 - sigx5_2 * sigx3_2)
  dw1pmu <- (0.5 * (dwv1) - 0.5 * (epsi_vu * pmusig))
  dw2pmu <- (dwv2 - epsi_vu2 * pmusig2)
  sigx7 <- (dw1pmu * sigx3_1 - dw2pmu * sigx3_2)
  Xsig1 <- sweep(2 * Xepsi_i, MARGIN = 1, STATS = pmusig2/ewv, FUN = "*") - sweep(X_iM,
    MARGIN = 1, STATS = 2 * (S * epsi_vu2/TT) * pmusig2/ewv, FUN = "*")
  Xsig2 <- sweep(X_iM, MARGIN = 1, STATS = S * dmusig2/ssq_v, FUN = "*")
  Xsig3 <- sweep(2 * Xepsi_i, MARGIN = 1, STATS = pmusig/ewv, FUN = "*") - sweep(X_iM,
    MARGIN = 1, STATS = 2 * (S * epsi_vu/TT) * pmusig/ewv, FUN = "*")
  Xsig4 <- sweep(X_iM, MARGIN = 1, STATS = S * dmusig/ssq_v, FUN = "*")
  Xsig5 <- sweep((0.5 * Xsig1 + Xsig2), MARGIN = 1, STATS = sigx3_2, FUN = "*") -
    sweep((0.5 * Xsig3 + Xsig4), MARGIN = 1, STATS = sigx3_1, FUN = "*")
  wusq <- (2/ewu_h - TT * epsi_vu2/ssq_v^2)
  sigx8 <- (1/ewu_h - TT * epsi_vu/ssq_v^2)
  sigx9 <- (ewv/(ewu_h * ssq_v^2))
  sigx10 <- sigx2_1 * epsi_vu/TT
  sigx11 <- dmusig * epsi_vu/ssq_v
  wutsq <- (wutt/(wutt)^2)
  Xsig6 <- sweep(X_iM, MARGIN = 1, STATS = (S * (1/(wutt) - epsi_vu/ssq_v^2)),
    FUN = "*")
  Xsig7 <- sweep(X_iM, MARGIN = 1, STATS = (S * (2/(wutt) - epsi_vu2/ssq_v^2)),
    FUN = "*")
  Xsig8 <- sweep(X_iM, MARGIN = 1, STATS = S * pmusig2, FUN = "*")
  Xsig9 <- sweep(X_iM, MARGIN = 1, STATS = S * (0.5 * (sigx11) + 0.5 * (pmusig -
    sigx11)), FUN = "*")
  Xsig10 <- sweep(X_iM, MARGIN = 1, STATS = S * epsivusq2 * dmusig2/ssq_v, FUN = "*")
  Xsig11 <- sweep(X_iM, MARGIN = 1, STATS = S * epsivusq1 * dmusig/ssq_v, FUN = "*")
  Xsig12 <- sweep(X_iM, MARGIN = 1, STATS = S * (sigx2_1 * epsi_vu/(TT * ewv) +
    0.5/ssq_v^2) * dmusig * ewv/sqrt(wvtt), FUN = "*")
  Xsig13 <- sweep(X_iM, MARGIN = 1, STATS = S * (epsi_vu2 * wusig/(TT * ewv) +
    0.5/ssq_v^2) * dmusig2 * ewv/sqrt(wvtt), FUN = "*")
  Xsig14 <- sweep(X_iM, MARGIN = 1, STATS = (S * epsi_vu2/TT), FUN = "*")
  Xsig15 <- sweep(X_iM, MARGIN = 1, STATS = (S * epsi_vu/TT), FUN = "*")
  Xsig16 <- sweep(X_iM, MARGIN = 1, STATS = 1/ssq_v, FUN = "*")
  Xsig17 <- sweep(X_iM, MARGIN = 1, STATS = S * pmusig, FUN = "*")
  Xsig18 <- sweep((2 * Xepsi_i - 2 * Xsig14), MARGIN = 1, STATS = S * dmusig2,
    FUN = "*")
  Xsig19 <- sweep((2 * Xepsi_i - 2 * Xsig14), MARGIN = 1, STATS = (dw2pmu/ewv),
    FUN = "*")
  Xsig20 <- sweep((2 * Xepsi_i - 2 * Xsig15), MARGIN = 1, STATS = (dw1pmu/ewv),
    FUN = "*")
  Xsig21 <- sweep(Xsig5, MARGIN = 1, STATS = TT * sigx7 * ewu_h/ttsq^2, FUN = "*")
  Xsig22 <- sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*")
  Xsig23 <- sweep((0.5 * Xsig19 + Xsig8), MARGIN = 1, STATS = sigx3_2, FUN = "*")
  Xsig24 <- sweep((0.5 * Xsig20 + Xsig9), MARGIN = 1, STATS = sigx3_1, FUN = "*")
  Xsig25 <- sweep((2 * Xepsi_i - 2 * Xsig15), MARGIN = 1, STATS = (sigx5_1/ewv),
    FUN = "*")
  Xsig26 <- sweep((2 * Xsig7 + 2 * Xsig22), MARGIN = 1, STATS = pmusig2, FUN = "*")
  Xsig27 <- sweep((2 * Xepsi_i - 2 * Xsig14), MARGIN = 1, STATS = (sigx5_2/ewv),
    FUN = "*")
  Xsig28 <- sweep((2 * Xsig6 + 2 * Xsig22), MARGIN = 1, STATS = pmusig, FUN = "*")
  Xsig29 <- sweep(Xsig5, MARGIN = 1, STATS = sigx6/sigx4, FUN = "*")
  Xsig30 <- sweep((0.5 * (Xsig28 - Xsig11) + Xsig12 - 0.5 * Xsig25), MARGIN = 1,
    STATS = sigx3_1, FUN = "*")
  Xsig31 <- sweep((0.5 * (Xsig26 - Xsig10) + Xsig13 - 0.5 * Xsig27), MARGIN = 1,
    STATS = sigx3_2, FUN = "*")
  Xsig32 <- sweep((Xsig30 - (Xsig29 + Xsig31)), MARGIN = 1, STATS = 1/sigx4, FUN = "*")
  Xsig33 <- sweep((Xsig23 - Xsig24), MARGIN = 1, STATS = 1/ttsq, FUN = "*")
  Xsig34 <- sweep((2 * Xepsi_i - 2 * Xsig15), MARGIN = 1, STATS = S * dmusig, FUN = "*")
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar + nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- 0.5 * (sapply(1:nXvar, function(x) crossprod(2 *
    Xsq[[x]], as.matrix(wHvar * pmusig2 * sigx3_2/ewv/sigx4))) - crossprod(sweep(X_iM,
    MARGIN = 1, STATS = wHvar * 2 * (S^2/TT) * pmusig2 * sigx3_2/ewv/sigx4, FUN = "*"),
    X_iM) - crossprod(sweep(Xsig18, MARGIN = 1, STATS = wHvar * sigx3_2/ewv/sigx4,
    FUN = "*"), Xsig16)) - (0.5 * crossprod(sweep((0.5 * Xsig1 + Xsig2), MARGIN = 1,
    STATS = wHvar * sigx3_2/ewv/sigx4, FUN = "*"), (2 * Xepsi_i - 2 * Xsig14)) +
    crossprod(sweep(X_iM, MARGIN = 1, STATS = wHvar * S^2 * epsi_vu2 * dmusig2/(TT^2 *
      sqrt(wvtt)) * sigx3_2/ewv/sigx4, FUN = "*"), X_iM)) - (0.5 * (sapply(1:nXvar,
    function(x) crossprod(2 * Xsq[[x]], as.matrix(wHvar * pmusig * sigx3_1/ewv/sigx4))) -
    crossprod(sweep(X_iM, MARGIN = 1, STATS = wHvar * 2 * (S^2/TT) * pmusig *
      sigx3_1/ewv/sigx4, FUN = "*"), X_iM) - crossprod(sweep(Xsig34, MARGIN = 1,
    STATS = wHvar * sigx3_1/ewv/sigx4, FUN = "*"), Xsig16)) - (0.5 * crossprod(sweep((0.5 *
    Xsig3 + Xsig4), MARGIN = 1, STATS = wHvar * sigx3_1/ewv/sigx4, FUN = "*"),
    (2 * Xepsi_i - 2 * Xsig15)) + crossprod(sweep(X_iM, MARGIN = 1, STATS = wHvar *
    S^2 * epsi_vu * dmusig/(TT^2 * sqrt(wvtt)) * sigx3_1/ewv/sigx4, FUN = "*"),
    X_iM))) - crossprod(sweep(Xsig5, MARGIN = 1, STATS = wHvar/sigx4^2, FUN = "*"),
    Xsig5)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod((Xsig33 - Xsig21),
    sweep(uHvar, MARGIN = 1, STATS = wHvar, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(Xsig32,
    sweep(vHvar, MARGIN = 1, STATS = wHvar, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.25 * (sigx11) - 0.5 * (0.5 * (sigx11) -
      0.5 * pmusig)) * ewv - 0.5 * (dw1pmu * epsi_vu/TT)) * sigx3_1 - sigx3_2 *
      (ewv * pmusig2 - epsi_vu2 * dw2pmu/TT))/(TT * sigx4 * ewu_h^2) - sigx7 *
      (dw1pmu * sigx3_1 + 0.5 * ttsq - dw2pmu * sigx3_2)/ttsq^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((0.5 * (0.5 *
    (epsivusq1 * dmusig * ewv/(wutt * sqrt(wvtt))) + 2 * (((0.25 * sigx9 - 0.5 *
    wutsq) * epsi_vu - 0.5 * (sigx2_1 * ewv/ewu_h)) * pmusig)) - (((0.25 * (ewv/ssq_v^2) +
    0.5 * (sigx10))/ewu_h - 0.5 * wutsq) * dwv1 + 0.5 * (sigx5_1 * epsi_vu/(wutt)))) *
    sigx3_1 - (sigx6 * sigx7/ttsq + (0.5 * (epsivusq2 * dmusig2 * ewv/(wutt *
    sqrt(wvtt)) + 2 * (((0.5 * sigx9 - wutt/(wutt)^2) * epsi_vu2 - wusig * ewv/ewu_h) *
    pmusig2)) - (((epsi_vu2 * wusig/TT + 0.5 * (ewv/ssq_v^2))/ewu_h - wutt/(wutt)^2) *
    dwv2 + sigx5_2 * epsi_vu2/(wutt))) * sigx3_2))/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (sigx5_1 * epsivusq1) + 0.5 * ((2 * ((sigx2_1/ewu_h - 0.5 * (sigx8 *
      epsi_vu/ssq_v^2)) * ewv) - epsilon_isq/ewv) * pmusig - sigx2_1 * epsivusq1 *
      dwv1) - (1/(wutt) - ((sigx2_1^2 + 0.5/ssq_v^2) * epsi_vu + 0.5 * (sigx8 *
      ewv/ssq_v^2) + 0.5 * sigx2_1)) * dwv1) * sigx3_1 - (sigx6^2/sigx4 + (0.5 *
      (sigx5_2 * epsivusq2) + 0.5 * ((2 * ((2 * (wusig/ewu_h) - 0.5 * (epsi_vu2 *
      wusq/ssq_v^2)) * ewv) - epsilon_isq/ewv) * pmusig2 - epsivusq2 * wusig *
      dwv2) - (2/(wutt) - ((wusig^2 + 0.5/ssq_v^2) * epsi_vu2 + 0.5 * (wusq *
      ewv/ssq_v^2) + 0.5 * wusig)) * dwv2) * sigx3_2))/sigx4, FUN = "*"), vHvar)
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
genexponormAlgOpt_pl81 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstgenexponorm_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    Initgenexpo <- start_st$initgenexpo
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(pgenexponormlike_pl81(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    -sum(pgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pgenexponormlike_pl81,
    grad = pgradgenexponormlike_pl81, hess = phessgenexponormlike_pl81, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessgenexponormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradgenexponormlike_pl81(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phessgenexponormlike_pl81(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessgenexponormlike_pl81(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- pgenexponormlike_pl81(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradgenexponormlike_pl81(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, Initgenexpo = Initgenexpo))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for genexpo-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
pgenexponormeff_pl81 <- function(object, level) {
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
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  mustar1 <- -(exp(Wv)/(TT * exp(Wu/2)) + object$S * epsilon_i/TT)
  mustar2 <- -(2 * exp(Wv)/(TT * exp(Wu/2)) + object$S * epsilon_i/TT)
  sigmastar <- sqrt(exp(Wv)/TT)
  a <- mustar1/sigmastar
  b <- mustar2/sigmastar
  A <- -1/2 * (epsilon_isq/exp(Wv) - a^2)
  B <- -1/2 * (epsilon_isq/exp(Wv) - b^2)
  u <- (exp(A) * (dnorm(a) * sigmastar + mustar1 * pnorm(a)) - exp(B) * (dnorm(b) *
    sigmastar + mustar2 * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- (exp(A) * exp(1/2 * sigmastar^2 - mustar1) * pnorm(a - sigmastar) -
      exp(B) * exp(1/2 * sigmastar^2 - mustar2) * pnorm(b - sigmastar))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal <- (exp(A) * exp(1/2 * sigmastar^2 + mustar1) * pnorm(a +
      sigmastar) - exp(B) * exp(1/2 * sigmastar^2 + mustar2) * pnorm(b + sigmastar))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    res <- data.frame(levels(pindex[, 1]), u = u, teJLMS = teJLMS, teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(levels(pindex[, 1]), u = u)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
