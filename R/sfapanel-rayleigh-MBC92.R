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
praynormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -exp(Wu) * S * giepsi/(exp(Wv) + gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  ll <- 1/2 * (mustar/sigmastar)^2 - epsilon_isq/(2 * exp(Wv)) - (TT - 1)/2 * log(2 *
    pi) - (TT - 1)/2 * Wv - 1/2 * Wu - 1/2 * log(exp(Wv) + gisq * exp(Wu)) +
    log(sigmastar * dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
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
pstraynorm_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, printInfo, tol, whichStart, initIter, initAlg) {
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
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "eta1", "eta2")
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
pgradraynormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
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
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  musig <- (S * ewu * giepsi/(sigmasq * sigmastar))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig, 0, 1)
  wvsq <- (1 - ewv/sigmasq)
  sigx1 <- (sigmastar/ewv - ewu/(sigmasq * sigmastar))
  sigx2 <- (pmusig + S * dmusig * sigx1 * giepsi)
  sigx3 <- (dmusig * sigmastar - S * ewu * pmusig * giepsi/sigmasq)
  sigx4 <- (S * giepsi/ewv - sigx2/sigx3)
  sigx5 <- (0.5 * (dmusig * ewv/sigmastar) - S * pmusig * giepsi)
  sigx6 <- ((1 - ewu * gisq/sigmasq) * ewv/sigmastar)
  sigx7 <- (1/(sigmasq * sigmastar) - (0.5 * sigx6 + sigmastar * gisq) * ewu/ssq)
  sigx8 <- (sigx5 * (1 - ewu * gisq/sigmasq)/sigx3 + S^2 * sigx7 * ewu * giepsi^2/sigmastar -
    0.5 * gisq)
  sigx9 <- (wvsq * dmusig/sigmastar)
  sigx10 <- (0.5 * sigx9 + S * pmusig * giepsi/sigmasq)
  sigx11 <- (0.5 * (wvsq * ewu/sigmastar) + sigmastar)
  sigx12 <- ((sigx10/sigx3 - S^2 * sigx11 * ewu * giepsi^2/(ssq * sigmastar)) *
    ewu - 0.5)
  sigx13 <- (sigmastar - 0.5 * (ewu * ewv/(sigmasq * sigmastar)))
  sigx14 <- (Zisq_epsi/(sigmasq * sigmastar) - ewu * sigx13 * giepsi * giZisq/ssq)
  sigx15 <- (0.5 * (dmusig * ewv * giZisq/(sigmasq * sigmastar)) + S^2 * dmusig *
    giepsi * sigx14)
  sigx16 <- (pmusig * giZisq/sigmasq + S * dmusig * sigx14)
  sigx17 <- (sigx15 * ewu + S * (pmusig * Zisq_epsi - ewu * sigx16 * giepsi))
  sigx18 <- (dmusig * ewv * giZi/(sigmasq * sigmastar))
  sigx19 <- (Zi_epsi/(sigmasq * sigmastar) - ewu * sigx13 * giepsi * giZi/ssq)
  sigx20 <- (pmusig * Zi_epsi - ewu * (pmusig * giZi/sigmasq + S * dmusig * sigx19) *
    giepsi)
  sigx21 <- (((0.5 * sigx18 + S^2 * dmusig * giepsi * sigx19) * ewu + S * sigx20)/sigx3 +
    0.5 * giZi)
  sigx22 <- (sigx17/sigx3 + 0.5 * giZisq)
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * ewu * sigx4/sigmasq, FUN = "*") -
    sweep(Xepsi_i, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = (sigx8 * ewu/sigmasq - 0.5), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = ((sigx12/sigmasq + 2 * (epsilon_isq/(2 * ewv)^2)) * ewv -
      0.5 * (TT - 1)), FUN = "*"), ewu * (S^2 * ewu * giepsi * sigx19/sigmastar -
    sigx21)/sigmasq, ewu * (S^2 * ewu * giepsi * sigx14/sigmastar - sigx22)/sigmasq)
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
phessraynormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
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
  Zisq <- as.numeric(tapply(2 * Zit^2, pindex[, 1], sum))
  Zicub <- as.numeric(tapply(2 * Zit^3, pindex[, 1], sum))
  Zifour <- as.numeric(tapply(2 * Zit^4, pindex[, 1], sum))
  XitZit <- sweep(-Xvar, MARGIN = 1, STATS = Zit, FUN = "*")
  XiZi <- apply(XitZit, 2, function(x) tapply(x, pindex[, 1], sum))
  XitZitsq <- sweep(-Xvar, MARGIN = 1, STATS = Zit^2, FUN = "*")
  XiZisq <- apply(XitZitsq, 2, function(x) tapply(x, pindex[, 1], sum))
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  musig <- (S * ewu * giepsi/(sigmasq * sigmastar))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig, 0, 1)
  wvsq <- (1 - ewv/sigmasq)
  ewusq <- (1 - ewu * gisq/sigmasq)
  sigx1 <- (sigmastar/ewv - ewu/(sigmasq * sigmastar))
  sigx2 <- (pmusig + S * dmusig * sigx1 * giepsi)
  sigx3 <- (dmusig * sigmastar - S * ewu * pmusig * giepsi/sigmasq)
  sigx4 <- (S * giepsi/ewv - sigx2/sigx3)
  sigx5 <- (0.5 * (dmusig * ewv/sigmastar) - S * pmusig * giepsi)
  sigx6 <- ((1 - ewu * gisq/sigmasq) * ewv/sigmastar)
  sigx7 <- (1/(sigmasq * sigmastar) - (0.5 * sigx6 + sigmastar * gisq) * ewu/ssq)
  sigx8 <- (sigx5 * (1 - ewu * gisq/sigmasq)/sigx3 + S^2 * sigx7 * ewu * giepsi^2/sigmastar -
    0.5 * gisq)
  sigx9 <- (wvsq * dmusig/sigmastar)
  sigx10 <- (0.5 * sigx9 + S * pmusig * giepsi/sigmasq)
  sigx11 <- (0.5 * (wvsq * ewu/sigmastar) + sigmastar)
  sigx12 <- ((sigx10/sigx3 - S^2 * sigx11 * ewu * giepsi^2/(ssq * sigmastar)) *
    ewu - 0.5)
  sigx13 <- (sigmastar - 0.5 * (ewu * ewv/(sigmasq * sigmastar)))
  sigx14 <- (Zisq_epsi/(sigmasq * sigmastar) - ewu * sigx13 * giepsi * giZisq/ssq)
  sigx15 <- (0.5 * (dmusig * ewv * giZisq/(sigmasq * sigmastar)) + S^2 * dmusig *
    giepsi * sigx14)
  sigx16 <- (pmusig * giZisq/sigmasq + S * dmusig * sigx14)
  sigx17 <- (sigx15 * ewu + S * (pmusig * Zisq_epsi - ewu * sigx16 * giepsi))
  sigx18 <- (dmusig * ewv * giZi/(sigmasq * sigmastar))
  sigx19 <- (Zi_epsi/(sigmasq * sigmastar) - ewu * sigx13 * giepsi * giZi/ssq)
  sigx20 <- (pmusig * Zi_epsi - ewu * (pmusig * giZi/sigmasq + S * dmusig * sigx19) *
    giepsi)
  sigx21 <- (((0.5 * sigx18 + S^2 * dmusig * giepsi * sigx19) * ewu + S * sigx20)/sigx3 +
    0.5 * giZi)
  sigx22 <- (sigx17/sigx3 + 0.5 * giZisq)
  sigx23 <- (S^2 * ewu * giepsi * sigx14/sigmastar - sigx22)
  sigx24 <- (0.5 * sigx6 + sigmastar * gisq)
  sigx25 <- ((0.5 * ewusq - 1.5)/(sigmasq * sigmastar) - 0.5 * sigx7)
  sigx26 <- ((sigx25 * ewv - 2 * (sigx24 * sigmasq * sigmastar * sigx13/ssq)) *
    ewu + 3 * sigmastar)
  sigx27 <- ((0.5 * ewusq - 0.5)/(sigmasq * sigmastar) - 0.5 * sigx7)
  sigx28 <- ((sigx27 * ewv - 2 * (sigx24 * sigmasq * sigmastar * sigx13/ssq)) *
    ewu + sigmastar)
  sigx29 <- (sigx28 * giepsi * giZi + sigx24 * Zi_epsi)
  sigx30 <- (sigx28 * giepsi * giZisq + sigx24 * Zisq_epsi)
  sigx31 <- (sigx30/ssq + S^2 * sigx7 * ewu * giepsi^2 * sigx14/(sigmasq * sigmastar))
  sigx32 <- (sigx29/ssq + S^2 * sigx7 * ewu * giepsi^2 * sigx19/(sigmasq * sigmastar))
  sigx33 <- (sigx24 * ewu/ssq - 1/(sigmasq * sigmastar))
  sigx34 <- (sigx24 * dmusig * ewv/ssq + S^2 * sigx7 * dmusig * giepsi^2/sigmasq)
  sigx35 <- ((dmusig + S^2 * sigx11 * dmusig * ewu^2 * ewv * giepsi^2/(ssq * sigmasq *
    sigmastar))/(sigmasq * sigmastar) - sigx11 * dmusig * ewv/ssq)
  sigx36 <- (S * sigx11 * dmusig * ewu * giepsi/ssq - pmusig/sigmasq)
  sigx37 <- ((0.5 * (wvsq/(sigmasq * sigmastar)) - 0.5 * (1/(sigmasq * sigmastar) -
    sigx11 * ewv/ssq)) * ewu - 2 * (sigx11 * sigmasq * sigmastar * sigx13/ssq))
  sigx38 <- ((0.5 * (sigx13/ssq) - 0.5/(sigmasq^2 * sigmastar)) * ewu^2 * ewv *
    giepsi * giZi + sigx13 * Zi_epsi)
  sigx39 <- (sigx38 * giZi + sigx13 * (giepsi * (Zisq - 2 * (ewu * sigmasq * sigmastar *
    sigx13 * giZi^2/ssq)) + giZi * Zi_epsi))
  sigx40 <- (sigx38 * giZisq + sigx13 * (giepsi * (Zicub - 2 * (ewu * sigmasq *
    sigmastar * sigx13 * giZi * giZisq/ssq)) + giZi * Zisq_epsi))
  sigx41 <- ((0.5 * (sigx13/ssq) - 0.5/(sigmasq^2 * sigmastar)) * ewu^2 * ewv *
    giepsi * giZisq + sigx13 * Zisq_epsi)
  sigx42 <- (sigx41 * giZisq + sigx13 * (giepsi * (Zifour - 2 * (ewu * sigmasq *
    sigmastar * sigx13 * giZisq^2/ssq)) + giZisq * Zisq_epsi))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 2, ncol = nXvar + nuZUvar +
    nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S^2 * (1/ewv - (((1 - S^2 * ewu * giepsi^2/(sigmasq * ewv)) * sigx1 - ewu/(sigmasq *
    sigmastar)) * dmusig + ewu * sigx2^2/(sigx3 * sigmasq))/sigx3) * ewu/sigmasq,
    FUN = "*"), Xgi) - sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar/(2 *
    ewv))))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * (((sigx5 * sigx2/sigx3 + 0.5 * (S * dmusig * giepsi/sigmastar)) *
      ewu/sigmasq - pmusig) * ewusq/sigx3 + 2 * (S * sigx7 * ewu * giepsi/sigmastar)) *
      ewu/sigmasq, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xepsi_i,
    MARGIN = 1, STATS = wHvar * 2/(2 * ewv)^2 * ewv, FUN = "*") + sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * (((sigx10 * sigx2/sigx3 - S * (0.5 * (wvsq/ewv) +
      1/sigmasq) * dmusig * giepsi/sigmastar) * ewu + pmusig)/(sigx3 * sigmasq) -
      2 * (S * sigx11 * ewu * giepsi/(ssq * sigmastar))) * ewu/sigmasq * ewv,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(XiZi, MARGIN = 1,
    STATS = wHvar/(sigmasq * sigmastar) * giepsi * S * ewu/sigmastar * S * ewu/sigmasq,
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * ewu * sigx13 * giZi/ssq *
    giepsi * S * ewu/sigmastar * S * ewu/sigmasq, FUN = "*") + sweep(Xgi, MARGIN = 1,
    STATS = wHvar * sigx19 * S * ewu/sigmastar * S * ewu/sigmasq, FUN = "*") -
    (sweep(Xgi, MARGIN = 1, STATS = wHvar * ((((0.5 * sigx18 + S^2 * dmusig *
      giepsi * sigx19) * ewu + S * sigx20) * sigx2/sigx3 - pmusig * giZi)/sigmasq +
      S * ((1 - S^2 * ewu * giepsi^2/(sigmasq * ewv)) * sigx19 + (ewu * (S^2 *
        giepsi * sigx19/ewv + giZi/(sigmasq * sigmastar)) * giepsi - (0.5 *
        (ewu * giepsi * giZi/sigmasq) + 2 * Zi_epsi)/sigmastar)/sigmasq +
        ewu * sigx13 * giepsi * giZi/ssq) * dmusig) * ewu/sigx3 * S * ewu/sigmasq,
      FUN = "*") + sweep(XiZi, MARGIN = 1, STATS = wHvar * pmusig/sigx3 * S *
      ewu/sigmasq, FUN = "*")))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(XiZisq, MARGIN = 1,
    STATS = wHvar/(sigmasq * sigmastar) * giepsi * S * ewu/sigmastar * S * ewu/sigmasq,
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * ewu * sigx13 * giZisq/ssq *
    giepsi * S * ewu/sigmastar * S * ewu/sigmasq, FUN = "*") + sweep(Xgi, MARGIN = 1,
    STATS = wHvar * sigx14 * S * ewu/sigmastar * S * ewu/sigmasq, FUN = "*") -
    (sweep(Xgi, MARGIN = 1, STATS = wHvar * ((sigx17 * sigx2/sigx3 - pmusig *
      giZisq)/sigmasq + S * ((1 - S^2 * ewu * giepsi^2/(sigmasq * ewv)) * sigx14 +
      (ewu * (S^2 * giepsi * sigx14/ewv + giZisq/(sigmasq * sigmastar)) * giepsi -
        (0.5 * (ewu * giepsi * giZisq/sigmasq) + 2 * Zisq_epsi)/sigmastar)/sigmasq +
      ewu * sigx13 * giepsi * giZisq/ssq) * dmusig) * ewu/sigx3 * S * ewu/sigmasq,
      FUN = "*") + sweep(XiZisq, MARGIN = 1, STATS = wHvar * pmusig/sigx3 *
      S * ewu/sigmasq, FUN = "*")))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (dmusig * ewv/sigmastar) + ewu * (S^2 *
      sigx7 * dmusig * giepsi^2 - (sigx5 * ewusq/sigx3 + gisq) * sigx5/sigmasq) -
      (0.5 * (0.5 * (ewusq * dmusig * ewv/sigmastar) + S^2 * sigx7 * dmusig *
        ewu * giepsi^2) + S * pmusig * giepsi)) * ewusq/sigx3 + ewu * (S^2 *
      (2/(sigmasq * sigmastar) - (((0.5 * (ewu * gisq/sigmasq) + 2 - 0.5 *
        (0.5 * ewusq + ewu * gisq/sigmasq)) * ewusq * ewv/sigmastar + (4 *
        gisq - 2 * (sigx24^2 * ewu * sigmasq/ssq)) * sigmastar) * ewu/ssq +
        0.5 * (ewusq * sigx7))) * giepsi^2/sigmastar - sigx8 * gisq/sigmasq) -
      0.5 * gisq) * ewu/sigmasq, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (((0.5 *
    (wvsq * dmusig) + 0.5 * ((dmusig * ewv * gisq/sigmasq - S^2 * wvsq * sigx7 *
    dmusig * ewu * giepsi^2/sigmastar) * ewu/sigmasq - 0.5 * (ewusq * wvsq *
    dmusig)))/sigmastar + (S * pmusig * giepsi - (sigx10 * sigx5 * ewusq/sigx3 +
    S * (pmusig * gisq/sigmasq + S * sigx7 * dmusig * giepsi) * giepsi) * ewu)/sigmasq)/sigx3 -
    (sigx12 * gisq/sigmasq + S^2 * (((0.5 * (ewusq * ewv/sigmasq) + 0.5 * ((ewu *
      gisq/sigmasq - 1) * ewv/sigmasq + 1 - 0.5 * (ewusq * wvsq)) + 0.5 * wvsq) *
      ewu/sigmastar + sigmastar)/(ssq * sigmastar) + sigx11 * (1/(ssq * sigmastar) -
      (0.5 * (ssq * ewusq/(sigmasq * sigmastar)) + 2 * (sigx24 * ewu)) * ewu *
        ewv/(ssq * sigmastar)^2)) * ewu * giepsi^2)) * ewu * ewv/sigmasq,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ewu * (ewu * (S^2 * (2 * (Zi_epsi/(sigmasq *
      sigmastar)) - ((sigx26 * giepsi * giZi + sigx24 * Zi_epsi) * ewu/ssq +
      0.5 * (ewusq * sigx19))) * giepsi/sigmastar - (((0.5 * (dmusig * ewv *
      giZi/sigmastar) - ((0.5 * sigx18 + S^2 * dmusig * giepsi * sigx19) *
      ewu + S * sigx20) * sigx5 * ewusq/sigx3)/sigmasq + S * (((pmusig * gisq/sigmasq +
      S * sigx7 * dmusig * giepsi) * giZi/sigmasq + S * sigx32 * dmusig) *
      ewu + S * sigx33 * dmusig * Zi_epsi - pmusig * giZi/sigmasq) * giepsi -
      (0.5 * (sigx34 * giZi) + S^2 * sigx32 * dmusig * giepsi) * ewu)/sigx3 +
      (S^2 * ewu * giepsi * sigx19/sigmastar - sigx21) * gisq/sigmasq)) - sigx21)/sigmasq,
    FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ewu * (ewu * (S^2 * (2 * (Zisq_epsi/(sigmasq *
      sigmastar)) - ((sigx26 * giepsi * giZisq + sigx24 * Zisq_epsi) * ewu/ssq +
      0.5 * (ewusq * sigx14))) * giepsi/sigmastar - (((0.5 * (dmusig * ewv *
      giZisq/sigmastar) - sigx17 * sigx5 * ewusq/sigx3)/sigmasq + S * (((pmusig *
      gisq/sigmasq + S * sigx7 * dmusig * giepsi) * giZisq/sigmasq + S * sigx31 *
      dmusig) * ewu + S * sigx33 * dmusig * Zisq_epsi - pmusig * giZisq/sigmasq) *
      giepsi - (0.5 * (sigx34 * giZisq) + S^2 * sigx31 * dmusig * giepsi) *
      ewu)/sigx3 + sigx23 * gisq/sigmasq)) - sigx22)/sigmasq, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((0.5 * (ewv * (S^2 * sigx11 * dmusig * ewu^2 * giepsi^2/(ssq * sigmastar) -
      dmusig)/sigmasq - 0.5 * (wvsq * dmusig)) + 0.5 * dmusig) * wvsq/sigmastar +
      (ewv * (S * sigx36 * giepsi - sigx10^2 * ewu/sigx3) + S * pmusig * giepsi)/sigmasq)/sigx3 -
      S^2 * (sigx11 * (1/(ssq * sigmastar) - (0.5 * (ssq * wvsq/(sigmasq *
        sigmastar)) + 2 * (sigx11 * ewv)) * ewu * ewv/(ssq * sigmastar)^2) +
        (0.5 * (ewv/sigmasq) - 0.5 * (0.5 * wvsq + ewv/sigmasq)) * wvsq *
          sigmasq/(ssq * ewv)) * ewu * giepsi^2) * ewu - (sigx12 * ewv/sigmasq +
      0.5))/sigmasq + (2 - 16 * (ewv^2/(2 * ewv)^2)) * epsilon_isq/(2 * ewv)^2) *
    ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar * ((((0.5 *
    (sigx35 * giZi) + S * (S * sigx11 * dmusig * Zi_epsi/ssq - sigx36 * giZi/sigmasq) *
    giepsi - ((0.5 * sigx18 + S^2 * dmusig * giepsi * sigx19) * ewu + S * sigx20) *
    sigx10/(sigx3 * sigmasq)) * ewv/sigx3 + S^2 * ((sigx37 * ewu * giepsi * giZi +
    sigx11 * Zi_epsi) * ewv/ssq + 0.5 * (wvsq * sigx19)) * giepsi/sigmastar) *
    ewu + ewv * (S^2 * ewu * giepsi * sigx19/sigmastar - sigx21)/sigmasq) * ewu/sigmasq),
    FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar * ((((0.5 *
    (sigx35 * giZisq) + S * (S * sigx11 * dmusig * Zisq_epsi/ssq - sigx36 * giZisq/sigmasq) *
    giepsi - sigx17 * sigx10/(sigx3 * sigmasq)) * ewv/sigx3 + S^2 * ((sigx37 *
    ewu * giepsi * giZisq + sigx11 * Zisq_epsi) * ewv/ssq + 0.5 * (wvsq * sigx14)) *
    giepsi/sigmastar) * ewu + ewv * sigx23/sigmasq) * ewu/sigmasq), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    (ewu * (ewu * (S^2 * ((0.5 * (giZi * sigx19/sigmasq) - sigx39/ssq) * ewu *
      giepsi + Zi_epsi * sigx19)/sigmastar - ((((0.5 * sigx18 + S^2 * dmusig *
      giepsi * sigx19) * ewu + S * sigx20)^2/(sigx3 * sigmasq) + 0.5 * (((dmusig *
      Zisq - S^2 * dmusig * ewu^2 * giepsi * giZi * sigx19/(sigmasq * sigmastar))/(sigmasq *
      sigmastar) - dmusig * ewu * sigx13 * giZi^2/ssq) * ewv) + S * (S * ((Zi_epsi -
      S^2 * ewu^2 * giepsi^2 * sigx19/(sigmasq * sigmastar)) * sigx19 - sigx39 *
      ewu * giepsi/ssq) * dmusig - (((pmusig * Zisq - ewu * (pmusig * giZi/sigmasq +
      S * dmusig * sigx19) * giZi)/sigmasq - S * (sigx39/ssq + S^2 * ewu *
      giepsi * sigx19^2/(sigmasq * sigmastar)) * dmusig * ewu) * giepsi + (2 *
      (S * dmusig * sigx19) + pmusig * giZi/sigmasq) * Zi_epsi)))/sigx3 + (S^2 *
      ewu * giepsi * sigx19/sigmastar - sigx21) * giZi/sigmasq)) - 0.5 * Zisq)/sigmasq))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    (ewu * (ewu * (S^2 * ((0.5 * (giZi * sigx14/sigmasq) - sigx40/ssq) * ewu *
      giepsi + Zi_epsi * sigx14)/sigmastar - ((((0.5 * sigx18 + S^2 * dmusig *
      giepsi * sigx19) * ewu + S * sigx20) * sigx17/(sigx3 * sigmasq) + 0.5 *
      (((dmusig * Zicub - S^2 * dmusig * ewu^2 * giepsi * giZisq * sigx19/(sigmasq *
        sigmastar))/(sigmasq * sigmastar) - dmusig * ewu * sigx13 * giZi *
        giZisq/ssq) * ewv) + S * (S * ((Zi_epsi - S^2 * ewu^2 * giepsi^2 *
      sigx19/(sigmasq * sigmastar)) * sigx14 - (sigx40 * ewu * giepsi/ssq +
      sigx19 * Zisq_epsi)) * dmusig - (((pmusig * Zicub - ewu * (pmusig * giZi/sigmasq +
      S * dmusig * sigx19) * giZisq)/sigmasq - S * (sigx40/ssq + S^2 * ewu *
      giepsi * sigx19 * sigx14/(sigmasq * sigmastar)) * dmusig * ewu) * giepsi +
      sigx16 * Zi_epsi)))/sigx3 + sigx23 * giZi/sigmasq)) - 0.5 * Zicub)/sigmasq))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    (ewu * (ewu * (S^2 * ((0.5 * (giZisq * sigx14/sigmasq) - sigx42/ssq) * ewu *
      giepsi + Zisq_epsi * sigx14)/sigmastar - ((sigx17^2/(sigx3 * sigmasq) +
      0.5 * (((dmusig * Zifour - S^2 * dmusig * ewu^2 * giepsi * giZisq * sigx14/(sigmasq *
        sigmastar))/(sigmasq * sigmastar) - dmusig * ewu * sigx13 * giZisq^2/ssq) *
        ewv) + S * (S * ((Zisq_epsi - S^2 * ewu^2 * giepsi^2 * sigx14/(sigmasq *
      sigmastar)) * sigx14 - sigx42 * ewu * giepsi/ssq) * dmusig - (((pmusig *
      Zifour - ewu * sigx16 * giZisq)/sigmasq - S * (sigx42/ssq + S^2 * ewu *
      giepsi * sigx14^2/(sigmasq * sigmastar)) * dmusig * ewu) * giepsi + (2 *
      (S * dmusig * sigx14) + pmusig * giZisq/sigmasq) * Zisq_epsi)))/sigx3 +
      sigx23 * giZisq/sigmasq)) - 0.5 * Zifour)/sigmasq))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for rayleigh-normal distribution
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
raynormAlgOpt_mbc92 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstraynorm_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    InitRay <- start_st$initRay
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(praynormlike_mbc92(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel Modified BC92 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(praynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = praynormlike_mbc92,
    grad = pgradraynormlike_mbc92, hess = phessraynormlike_mbc92, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(praynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(praynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(praynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(praynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessraynormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradraynormlike_mbc92(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phessraynormlike_mbc92(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessraynormlike_mbc92(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- praynormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradraynormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
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
praynormeff_mbc92 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  eta1 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  eta2 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 2]
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -exp(Wu) * object$S * giepsi/(exp(Wv) + gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  u <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 + sigmastar^2) *
    pnorm(mustar/sigmastar))/(sigmastar * dnorm(mustar/sigmastar) + mustar *
    pnorm(mustar/sigmastar))
  m <- ifelse(mustar/2 + sqrt(sigmastar^2 + mustar^2/4) > 0, mustar/2 + sqrt(sigmastar^2 +
    mustar^2/4), 0)
  res <- data.frame(levels(pindex[, 1]), u = u, m = m, mustar = mustar, sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  res$m <- res$m * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teMO <- exp(-res$m)
    res$teBC <- exp(-res$mustar * git + res$sigmastar^2 * git^2/2) * (res$sigmastar *
      git * dnorm(res$mustar/res$sigmastar - res$sigmastar * git) + (res$mustar *
      git - res$sigmastar^2 * git^2) * pnorm(res$mustar/res$sigmastar - res$sigmastar *
      git))/(res$sigmastar * git * dnorm(res$mustar/res$sigmastar) + res$mustar *
      git * pnorm(res$mustar/res$sigmastar))
    res$teBC_reciprocal <- exp(res$mustar * git + res$sigmastar^2 * git^2/2) *
      (res$sigmastar * git * dnorm(res$mustar/res$sigmastar + res$sigmastar *
        git) + (res$mustar * git + res$sigmastar^2 * git^2) * pnorm(res$mustar/res$sigmastar +
        res$sigmastar * git))/(res$sigmastar * git * dnorm(res$mustar/res$sigmastar) +
      res$mustar * git * pnorm(res$mustar/res$sigmastar))
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
