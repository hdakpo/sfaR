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
phalfnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  ll <- log(2) - TT/2 * log(2 * pi) + pnorm(mustar/sigmastar, log.p = TRUE) + log(sigmastar) -
    1/2 * (epsilon_isq/exp(Wv) - (mustar/sigmastar)^2) - TT/2 * log(exp(Wv)) -
    1/2 * log(exp(Wu))
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
psthalfnorm_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA + halfnormal-normal distribution...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "eta1", "eta2")
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
pgradhalfnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar) {
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
  sigx1 <- (S * ewu * giepsi/(sigmasq * sigmastar) - dmusig/pmusig)
  wvsq <- (1 - ewv/sigmasq)
  sigx2 <- (wvsq * ewu/sigmastar)
  sigx3 <- (sigmastar - 0.5 * (ewu * ewv/(sigmasq * sigmastar)))
  wu2sq <- (1 - ewu * gisq/sigmasq)
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (1/(sigmasq * sigmastar) - sigx4 * ewu/ssq)
  sps <- (sigmasq * pmusig * sigmastar)
  sqspmu <- (ssq * pmusig)
  sigx6 <- (ssq * sigmasq * sigmastar)
  sigx7 <- (Zi_epsi/(sigmasq * sigmastar) - ewu * sigx3 * giepsi * giZi/ssq)
  sigx8 <- (Zisq_epsi/(sigmasq * sigmastar) - ewu * sigx3 * giepsi * giZisq/ssq)
  sigx9 <- (0.5 * sigx2 + sigmastar)
  gradll <- cbind(-(sweep(Xepsi_i, MARGIN = 1, STATS = 0.5/ewv, FUN = "*") - sweep(Xgi,
    MARGIN = 1, STATS = (S^2 * ewu * giepsi/sigmasq)/ewv, FUN = "*") + sweep(Xgi,
    MARGIN = 1, STATS = S * dmusig * ewu/sps, FUN = "*")), sweep(uHvar, MARGIN = 1,
    STATS = (0.5 * wu2sq + S * sigx5 * ewu * sigx1 * giepsi - 0.5), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (0.5 * wvsq + S * sigx9 * dmusig * ewu *
      ewv * giepsi/sqspmu - (0.5 * (2 * (S^2 * sigx9 * ewu^2 * ewv * giepsi^2/sigx6) -
      epsilon_isq/ewv) + 0.5 * TT)), FUN = "*"), ewu * (S * sigx1 * sigx7 -
      0.5 * (giZi/sigmasq)), ewu * (S * sigx1 * sigx8 - 0.5 * (giZisq/sigmasq)))
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
phesshalfnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar) {
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
  sps <- (sigmasq * pmusig * sigmastar)
  sqspmu <- (ssq * pmusig)
  ewvsqx2 <- (wvsq/(sigmasq * sigmastar))
  sigx1 <- (S * ewu * giepsi/(sigmasq * sigmastar) - dmusig/pmusig)
  sigx2 <- (wvsq * ewu/sigmastar)
  sigx3 <- (sigmastar - 0.5 * (ewu * ewv/(sigmasq * sigmastar)))
  wu2sq <- (1 - ewu * gisq/sigmasq)
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (1/(sigmasq * sigmastar) - sigx4 * ewu/ssq)
  sigx6 <- (ssq * sigmasq * sigmastar)
  sigx7 <- (Zi_epsi/(sigmasq * sigmastar) - ewu * sigx3 * giepsi * giZi/ssq)
  sigx8 <- (Zisq_epsi/(sigmasq * sigmastar) - ewu * sigx3 * giepsi * giZisq/ssq)
  sigx9 <- (0.5 * sigx2 + sigmastar)
  sigx10 <- (0.5 + 0.5 * (ewv/sigmasq) - 0.5 * (0.5 * wvsq + ewv/sigmasq))
  sigx11 <- (ewu * gisq/sigmasq)
  sigx12 <- ((ewu * gisq/sigmasq - 1) * ewv/sigmasq + 1 - 0.5 * (wu2sq * wvsq))
  sigx13 <- (1/(sigmasq * sigmastar) - sigx9 * ewv/ssq)
  sigx14 <- (0.5 * ewvsqx2 - 0.5 * sigx13) * ewu
  sigx15 <- (sigx9 * sigmasq * sigmastar * sigx3/ssq)
  sigx16 <- (dmusig/pmusig - S * ewu * giepsi/(sigmasq * sigmastar))
  sigx17 <- (1 - dmusig * sigx16/pmusig)
  sigx18 <- (1 + dmusig * sigx1/pmusig)
  sigx19 <- (1/sigmastar - dmusig * (dmusig/(pmusig * sigmastar) - S * giepsi/ewv)/pmusig)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 2, ncol = nXvar + nuZUvar +
    nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- -(0.5 * (sapply(1:nXvar, function(x) crossprod(Xsq[[x]],
    as.matrix(wHvar/ewv))) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    2 * (S^2 * ewu/sigmasq)/ewv, FUN = "*"), Xgi)) + crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S^2 * dmusig * (dmusig/sps^2 - S * giepsi/(sigmasq^2 * ewv *
      pmusig * sigmastar)) * ewu^2, FUN = "*"), Xgi))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * sigx5 * ewu * (S * (2/sigmastar - dmusig * (dmusig/(pmusig *
      sigmastar) - S * giepsi/ewv)/pmusig) * ewu * giepsi/sigmasq - dmusig/pmusig),
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * ((ewv - S^2 * ewu * giepsi^2/sigmasq)/sqspmu +
      S * ssq * dmusig * ewu * ewv * giepsi/(sqspmu^2 * sigmasq * sigmastar)) *
      sigx9 * dmusig * ewu, FUN = "*") - 0.5 * (sweep(Xgi, MARGIN = 1, STATS = wHvar *
    4 * (S^2 * sigx9 * ewu^2 * ewv * giepsi/sigx6), FUN = "*") - sweep(Xepsi_i,
    MARGIN = 1, STATS = wHvar/ewv, FUN = "*")), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(XiZi, MARGIN = 1,
    STATS = wHvar * sigx1/(sigmasq * sigmastar) * S * ewu, FUN = "*") - sweep(Xgi,
    MARGIN = 1, STATS = wHvar * sigx1 * ewu * sigx3 * giZi/ssq * S * ewu, FUN = "*") +
    sweep(Xgi, MARGIN = 1, STATS = wHvar * S * sigx19 * ewu * sigx7/sigmasq *
      S * ewu, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(XiZisq, MARGIN = 1,
    STATS = wHvar * sigx1/(sigmasq * sigmastar) * S * ewu, FUN = "*") - sweep(Xgi,
    MARGIN = 1, STATS = wHvar * sigx1 * ewu * sigx3 * giZisq/ssq * S * ewu, FUN = "*") +
    sweep(Xgi, MARGIN = 1, STATS = wHvar * S * sigx19 * ewu * sigx8/sigmasq *
      S * ewu, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ewu * (S * ((1/(sigmasq * sigmastar) - ((0.5 *
      sigx11 + 1.5 - 0.5 * (0.5 * wu2sq + ewu * gisq/sigmasq)) * wu2sq * ewv/sigmastar +
      (3 * gisq - 2 * (sigx4^2 * ewu * sigmasq/ssq)) * sigmastar) * ewu/ssq) *
      sigx1 + S * sigx17 * sigx5^2 * ewu * giepsi) * giepsi - 0.5 * (wu2sq *
      gisq/sigmasq)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (0.5 * (gisq/sigmasq^2) +
    S * (((((0.5 * (wu2sq * ewv) - S^2 * sigx9 * sigx5 * ewu * giepsi^2)/sigmasq +
      0.5 * sigx12 + 0.5 * wvsq) * ewu/sigmastar + sigmastar)/sqspmu - sigx9 *
      (2 * (sigx4 * sigmasq * pmusig * sigmastar) - S * ssq * sigx5 * dmusig *
        giepsi) * ewu/sqspmu^2) * dmusig - S * (((0.5 * (wu2sq * ewv/sigmasq) +
      0.5 * sigx12) * ewu/sigmastar + 2 * sigx9)/sigx6 - ((ssq * gisq + 2 *
      (sigx4 * sigmasq^2 * sigmastar)) * sigmastar + 0.5 * (ssq * wu2sq * ewv/sigmastar)) *
      sigx9 * ewu/sigx6^2) * ewu * giepsi) * giepsi) * ewu * ewv, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (gisq * giZi/sigmasq^2) + S * (S * sigx17 *
      sigx5 * giepsi * sigx7 - (((((0.5 * wu2sq - 0.5)/(sigmasq * sigmastar) -
      0.5 * sigx5) * ewv - 2 * (sigx4 * sigmasq * sigmastar * sigx3/ssq)) *
      ewu + sigmastar) * giepsi * giZi + sigx4 * Zi_epsi) * sigx1/ssq)) * ewu +
      S * sigx1 * sigx7 - 0.5 * (giZi/sigmasq)) * ewu, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (gisq * giZisq/sigmasq^2) + S * (S *
      sigx17 * sigx5 * giepsi * sigx8 - (((((0.5 * wu2sq - 0.5)/(sigmasq *
      sigmastar) - 0.5 * sigx5) * ewv - 2 * (sigx4 * sigmasq * sigmastar *
      sigx3/ssq)) * ewu + sigmastar) * giepsi * giZisq + sigx4 * Zisq_epsi) *
      sigx1/ssq)) * ewu + S * sigx1 * sigx8 - 0.5 * (giZisq/sigmasq)) * ewu,
    FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (ewv * (S * (((sigx10 * wvsq + S^2 * sigx9^2 * ewu * ewv * giepsi^2/(ssq *
      sigmasq)) * ewu/sigmastar + sigmastar)/sqspmu - sigx9^2 * (2 * sps +
      S * dmusig * ewu * giepsi) * ewv/sqspmu^2) * dmusig * ewu * giepsi -
      0.5 * (wvsq/sigmasq)) - 0.5 * (2 * (S^2 * ((sigx10 * wvsq * ewu/sigmastar +
      sigmastar)/sigx6 - ((ssq + 2 * (sigx9 * sigmasq^2 * sigmastar)) * sigmastar +
      0.5 * (ssq * wvsq * ewu/sigmastar)) * sigx9 * ewv/sigx6^2) * ewu^2 *
      ewv * giepsi^2) + epsilon_isq/ewv)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (0.5 *
    (giZi/sigmasq^2) - S * (((sigx14 - 2 * sigx15) * ewu * giepsi * giZi + sigx9 *
    Zi_epsi) * sigx1 + S * sigx9 * sigx18 * ewu * giepsi * sigx7)/ssq) * ewu *
    ewv, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (0.5 *
    (giZisq/sigmasq^2) - S * (((sigx14 - 2 * sigx15) * ewu * giepsi * giZisq +
    sigx9 * Zisq_epsi) * sigx1 + S * sigx9 * sigx18 * ewu * giepsi * sigx8)/ssq) *
    ewu * ewv, FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    (ewu * (S * ewu * (S * sigx17 * sigx7^2 - (((0.5 * (sigx3/ssq) - 0.5/(sigmasq^2 *
      sigmastar)) * ewu^2 * ewv * giepsi * giZi + sigx3 * Zi_epsi) * giZi +
      sigx3 * (giepsi * (Zisq - 2 * (ewu * sigmasq * sigmastar * sigx3 * giZi^2/ssq)) +
        giZi * Zi_epsi)) * sigx1/ssq) - 0.5 * ((Zisq - ewu * giZi^2/sigmasq)/sigmasq))))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    (ewu * (S * ewu * (S * sigx17 * sigx7 * sigx8 - (((0.5 * (sigx3/ssq) - 0.5/(sigmasq^2 *
      sigmastar)) * ewu^2 * ewv * giepsi * giZi + sigx3 * Zi_epsi) * giZisq +
      sigx3 * (giepsi * (Zicub - 2 * (ewu * sigmasq * sigmastar * sigx3 * giZi *
        giZisq/ssq)) + giZi * Zisq_epsi)) * sigx1/ssq) - 0.5 * ((Zicub -
      ewu * giZi * giZisq/sigmasq)/sigmasq))))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    (ewu * (S * ewu * (S * sigx17 * sigx8^2 - (((0.5 * (sigx3/ssq) - 0.5/(sigmasq^2 *
      sigmastar)) * ewu^2 * ewv * giepsi * giZisq + sigx3 * Zisq_epsi) * giZisq +
      sigx3 * (giepsi * (Zifour - 2 * (ewu * sigmasq * sigmastar * sigx3 *
        giZisq^2/ssq)) + giZisq * Zisq_epsi)) * sigx1/ssq) - 0.5 * ((Zifour -
      ewu * giZisq^2/sigmasq)/sigmasq))))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
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
halfnormAlgOpt_mbc92 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psthalfnorm_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    InitHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(phalfnormlike_mbc92(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    -sum(phalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradhalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = phalfnormlike_mbc92,
    grad = pgradhalfnormlike_mbc92, hess = phesshalfnormlike_mbc92, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(phalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(phalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesshalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(phalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradhalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesshalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(phalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradhalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phesshalfnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradhalfnormlike_mbc92(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phesshalfnormlike_mbc92(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesshalfnormlike_mbc92(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- phalfnormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradhalfnormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, InitHalf = InitHalf))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for halfnormal-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
phalfnormeff_mbc92 <- function(object, level) {
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
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) * sigmastar
  m <- ifelse(mustar > 0, mustar, 0)
  res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB, uUB = uUB, m = m, mustar = mustar,
    sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  res$m <- res$m * git
  res$uLB <- res$uLB * git
  res$uUB <- res$uUB * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teMO <- exp(-res$m)
    res$teBC <- exp(1/2 * res$sigmastar^2 * git^2 - res$mustar * git) * pnorm(res$mustar/res$sigmastar -
      res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
    res$teBCLB <- exp(-res$uUB)
    res$teBCUB <- exp(-res$uLB)
    res$teBC_reciprocal <- exp(1/2 * res$sigmastar^2 * git^2 + res$mustar * git) *
      pnorm(res$mustar/res$sigmastar + res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
