################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Two types: - Battese and Coelli 1992 - alternative (mbc92)                   #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = 1 + eta1 * (t - T) + eta2 * (t - T)^2                  #
# Convolution: exponential - normal                                            #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for exponential-normal distribution
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
pexponormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x -
    1))) + eta2 * (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -(exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  ll <- 1/2 * Wu - (TT - 1)/2 * Wv - (TT - 1)/2 * log(2 * pi) -
    1/2 * log(gisq) - epsilon_isq/(2 * exp(Wv)) + 1/2 * (mustar/sigmastar)^2 +
    pnorm(mustar/sigmastar, log.p = TRUE)
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for exponential-normal distribution
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
pstexponorm_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, printInfo, tol,
  whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstexponorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initExpo <- NULL
  } else {
    cat("Initialization: SFA + exponential-normal distribution...\n")
    initExpo <- maxLik::maxLik(logLik = cexponormlike, start = cstexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradexponormlike,
      method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initExpo$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "eta1",
    "eta2")
  return(list(StartVal = StartVal, initExpo = initExpo))
}

# Gradient of the likelihood function ----------
#' gradient for exponential-normal distribution
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
pgradexponormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x -
    1))) + eta2 * (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zit <- unlist(lapply(TT, FUN = function(x) rev(-(0:(x - 1)))))
  gitZit <- 2 * (Zit * (1 + Zit * (eta1 + eta2 * Zit)))
  giZi <- as.numeric(tapply(gitZit, pindex[, 1], sum))
  gitZitsq <- 2 * (Zit^2 * (1 + Zit * (eta1 + eta2 * Zit)))
  giZisq <- as.numeric(tapply(gitZitsq, pindex[, 1], sum))
  Zit_epsit <- Zit * epsilon_it
  Zi_epsi <- as.numeric(tapply(Zit_epsit, pindex[, 1], sum))
  Zitsq_epsit <- Zit^2 * epsilon_it
  Zisq_epsi <- as.numeric(tapply(Zitsq_epsit, pindex[, 1],
    sum))

  ewv <- exp(Wv)
  ewu_h <- exp(Wu/2)
  sigmastar <- sqrt(ewv/gisq)
  gisig <- (sigmastar * gisq)
  euvg <- (ewv/ewu_h + S * giepsi)
  dmusig <- dnorm(-(euvg/gisig), 0, 1)
  pmusig <- pnorm(-(euvg/gisig))
  sigx1 <- (euvg/ewv - dmusig/(pmusig * sigmastar))
  sigx2 <- (dmusig * ewv/(pmusig * sigmastar))
  sigx3 <- (1/ewu_h - dmusig/(pmusig * sigmastar))
  sigx4 <- (1/(ewu_h * gisq) - 0.5 * (euvg/gisig^2))
  sigx5 <- (sigx3 * ewv + S * giepsi)
  sigx6 <- (euvg/gisig - dmusig/pmusig)
  sigx7 <- (sigmastar - 0.5 * (ewv/gisig))
  sigx8 <- (S * Zi_epsi/gisig - euvg * sigx7 * giZi/gisig^2)
  sigx9 <- (S * Zisq_epsi/gisig - euvg * sigx7 * giZisq/gisig^2)
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * sigx1/gisq,
    FUN = "*") - sweep(Xepsi_i, MARGIN = 1, STATS = 1/(2 *
    ewv), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ((0.5 *
    sigx2 - 0.5 * euvg)/(ewu_h * gisq) + 0.5), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (sigx5 * sigx4 + 2 *
      (ewv * epsilon_isq/(2 * ewv)^2) - 0.5 * (TT - 1)),
      FUN = "*"), sigx6 * sigx8 - 0.5 * (giZi/gisq), sigx6 *
      sigx9 - 0.5 * (giZisq/gisq))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for exponential-normal distribution
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
phessexponormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x -
    1))) + eta2 * (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zit <- unlist(lapply(TT, FUN = function(x) rev(-(0:(x - 1)))))
  gitZit <- 2 * (Zit * (1 + Zit * (eta1 + eta2 * Zit)))
  giZi <- as.numeric(tapply(gitZit, pindex[, 1], sum))
  gitZitsq <- 2 * (Zit^2 * (1 + Zit * (eta1 + eta2 * Zit)))
  giZisq <- as.numeric(tapply(gitZitsq, pindex[, 1], sum))
  Zit_epsit <- Zit * epsilon_it
  Zi_epsi <- as.numeric(tapply(Zit_epsit, pindex[, 1], sum))
  Zitsq_epsit <- Zit^2 * epsilon_it
  Zisq_epsi <- as.numeric(tapply(Zitsq_epsit, pindex[, 1],
    sum))
  Zisq <- as.numeric(tapply(2 * Zit^2, pindex[, 1], sum))
  Zicub <- as.numeric(tapply(2 * Zit^3, pindex[, 1], sum))
  Zifour <- as.numeric(tapply(2 * Zit^4, pindex[, 1], sum))
  XitZit <- sweep(-Xvar, MARGIN = 1, STATS = Zit, FUN = "*")
  XiZi <- apply(XitZit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  XitZitsq <- sweep(-Xvar, MARGIN = 1, STATS = Zit^2, FUN = "*")
  XiZisq <- apply(XitZitsq, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) tapply(x, pindex[,
      1], sum))
  }
  ewv <- exp(Wv)
  ewu_h <- exp(Wu/2)
  sigmastar <- sqrt(ewv/gisq)
  gisig <- (sigmastar * gisq)
  euvg <- (ewv/ewu_h + S * giepsi)
  dmusig <- dnorm(-(euvg/gisig), 0, 1)
  pmusig <- pnorm(-(euvg/gisig))
  sigx1 <- (euvg/ewv - dmusig/(pmusig * sigmastar))
  sigx2 <- (dmusig * ewv/(pmusig * sigmastar))
  sigx3 <- (1/ewu_h - dmusig/(pmusig * sigmastar))
  sigx4 <- (1/(ewu_h * gisq) - 0.5 * (euvg/gisig^2))
  sigx5 <- (sigx3 * ewv + S * giepsi)
  sigx6 <- (euvg/gisig - dmusig/pmusig)
  sigx7 <- (sigmastar - 0.5 * (ewv/gisig))
  sigx8 <- (S * Zi_epsi/gisig - euvg * sigx7 * giZi/gisig^2)
  sigx9 <- (S * Zisq_epsi/gisig - euvg * sigx7 * giZisq/gisig^2)
  sigx10 <- ((0.5 * euvg - 0.5 * sigx2) * dmusig/pmusig + 0.5 *
    (ewv/sigmastar))
  sigx11 <- (0.5 * (euvg/(pmusig * sigmastar)) - 0.5 * (dmusig *
    ewv/(pmusig * sigmastar)^2))
  sigx12 <- (ewv/sigmastar - ((dmusig/(pmusig * sigmastar) -
    1/ewu_h) * ewv - S * giepsi) * dmusig/pmusig)
  sigx13 <- (((0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/gisig^2)))/sigmastar -
    sigx7 * gisq/gisig^2) * euvg + sigx7/ewu_h)
  sigx14 <- (0.5 * (sigx7/gisig^2) - 0.5/(sigmastar * gisq^2))
  sigx15 <- (sigx14 * ewv * euvg * giZi + S * sigx7 * Zi_epsi)
  sigx16 <- (sigx14 * ewv * euvg * giZisq + S * sigx7 * Zisq_epsi)
  sigx17 <- (dmusig/pmusig - euvg/gisig)
  sigx18 <- (1/sigmastar - dmusig * (dmusig/(pmusig * sigmastar) -
    euvg/ewv)/pmusig)
  sigx19 <- (dmusig/(pmusig * sigmastar)^2 - euvg/(ewv * pmusig *
    sigmastar))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 2, ncol = nXvar +
    nuZUvar + nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * (1/ewv - dmusig * sigx19/gisq)/gisq,
    FUN = "*"), Xgi) - sapply(1:nXvar, function(x) crossprod(Xsq[[x]],
    as.matrix(wHvar/(2 * ewv))))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * (0.5 * (dmusig * (dmusig *
      ewv/(pmusig * sigmastar)^2 - euvg/(pmusig * sigmastar))/gisq) -
      0.5)/(ewu_h * gisq), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xepsi_i, MARGIN = 1, STATS = wHvar *
    2 * (ewv/(2 * ewv)^2), FUN = "*") + sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * ((1 - dmusig * sigx19 * ewv/gisq) *
      sigx4 - 0.5 * (sigx5/gisig^2)), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- S * colSums(sweep(XiZi,
    MARGIN = 1, STATS = wHvar/gisig * sigx6, FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx7 * giZi/gisig^2 *
      sigx6, FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    sigx18 * sigx8/gisq, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2)] <- S * colSums(sweep(XiZisq,
    MARGIN = 1, STATS = wHvar/gisig * sigx6, FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx7 * giZisq/gisig^2 *
      sigx6, FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    sigx18 * sigx9/gisq, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((0.25 + 0.5 * (sigx11 * dmusig/gisq)) * ewv/(ewu_h^2 *
      gisq) - 0.5 * ((0.5 * sigx2 - 0.5 * euvg) * ewu_h *
      gisq/(ewu_h * gisq)^2)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx5 * (0.25 * (ewv/(ewu_h *
      gisig^2)) - 0.5 * (ewu_h * gisq/(ewu_h * gisq)^2)) -
      (sigx11 * dmusig/gisq + 0.5) * sigx4 * ewv/ewu_h),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (sigx6 * ewv * sigx7 * giZi/gisig^2) - sigx10 *
      sigx8/gisq)/ewu_h, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 2)] <- colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (sigx6 * ewv * sigx7 * giZisq/gisig^2) - sigx10 *
      sigx9/gisq)/ewu_h, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx4 * euvg - 1)/(pmusig *
      sigmastar) + (0.5 * (pmusig/gisig) - sigx4 * dmusig) *
      ewv/(pmusig * sigmastar)^2) * dmusig + 1/ewu_h) *
      sigx4 + 2 * ((1 - 8 * (ewv^2/(2 * ewv)^2)) * epsilon_isq/(2 *
      ewv)^2) - 0.5 * (sigx5 * (1/ewu_h - euvg * gisq/gisig^2)/gisig^2)) *
      ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx4 * sigx12 * sigx8 -
      (sigx13 * giZi + 0.5 * (S * Zi_epsi/sigmastar)) *
        sigx6 * ewv/gisig^2), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx4 * sigx12 * sigx9 -
      (sigx13 * giZisq + 0.5 * (S * Zisq_epsi/sigmastar)) *
        sigx6 * ewv/gisig^2), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar +
    nvZVvar + 1)] <- sum(wHvar * ((S * Zi_epsi/gisig - (euvg *
    sigx7 * giZi/gisig^2 + dmusig * sigx17 * sigx8/pmusig)) *
    sigx8 - ((sigx15 * giZi + (euvg * (Zisq - 2 * (sigmastar *
    sigx7 * gisq * giZi^2/gisig^2)) + S * giZi * Zi_epsi) *
    sigx7) * sigx6/gisig^2 + 0.5 * ((Zisq - giZi^2/gisq)/gisq))))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * ((S * Zi_epsi/gisig - (euvg *
    sigx7 * giZi/gisig^2 + dmusig * sigx17 * sigx8/pmusig)) *
    sigx9 - ((sigx15 * giZisq + (euvg * (Zicub - 2 * (sigmastar *
    sigx7 * gisq * giZi * giZisq/gisig^2)) + S * giZi * Zisq_epsi) *
    sigx7) * sigx6/gisig^2 + 0.5 * ((Zicub - giZi * giZisq/gisq)/gisq))))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * ((S * Zisq_epsi/gisig -
    (euvg * sigx7 * giZisq/gisig^2 + dmusig * sigx17 * sigx9/pmusig)) *
    sigx9 - ((sigx16 * giZisq + (euvg * (Zifour - 2 * (sigmastar *
    sigx7 * gisq * giZisq^2/gisig^2)) + S * giZisq * Zisq_epsi) *
    sigx7) * sigx6/gisig^2 + 0.5 * ((Zifour - giZisq^2/gisq)/gisq))))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for exponential-normal distribution
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
#' @param whichStart strategy to get starting values
#' @param itermax maximum iteration
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik 
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
exponormAlgOpt_mbc92 <- function(start, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar,
  Yvar, Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo,
  whichStart, initIter, initAlg, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstexponorm_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar_c, tol = tol, whichStart = whichStart,
      initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    InitExpo <- start_st$initExpo
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(pexponormlike_mbc92(startVal, nXvar = nXvar,
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
  cat("SFA Panel Modified BC92 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) {
      -sum(pexponormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradexponormlike_mbc92(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pexponormlike_mbc92,
    grad = pgradexponormlike_mbc92, hess = phessexponormlike_mbc92,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pexponormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradexponormlike_mbc92(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pexponormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradexponormlike_mbc92(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessexponormlike_mbc92(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) {
      -sum(pexponormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradexponormlike_mbc92(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessexponormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(pexponormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradexponormlike_mbc92(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phessexponormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradexponormlike_mbc92(mleObj$par,
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
      mleObj$hessian <- phessexponormlike_mbc92(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessexponormlike_mbc92(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- pexponormlike_mbc92(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradexponormlike_mbc92(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, InitExpo = InitExpo))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for exponential-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
pexponormeff_mbc92 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  eta1 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  eta2 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    2]
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x -
    1))) + eta2 * (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -(exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
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
    res$teBC <- exp(-res$mustar * git + 1/2 * res$sigmastar^2 *
      git^2) * pnorm(res$mustar/res$sigmastar - res$sigmastar *
      git)/pnorm(res$mustar/res$sigmastar)
    res$teBCLB <- exp(-res$uUB)
    res$teBCUB <- exp(-res$uLB)
    teBC_reciprocal <- exp(res$mustar * git + 1/2 * res$sigmastar^2 *
      git^2) * pnorm(res$mustar/res$sigmastar + res$sigmastar *
      git)/pnorm(res$mustar/res$sigmastar)
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
